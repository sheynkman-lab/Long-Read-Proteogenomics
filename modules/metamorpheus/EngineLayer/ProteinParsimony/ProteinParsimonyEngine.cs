﻿using EngineLayer.ProteinParsimony;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer
{
    public class ProteinParsimonyEngine : MetaMorpheusEngine
    {
        /// <summary>
        /// All peptides meeting the prefiltering criteria for parsimony (e.g., peptides from non-ambiguous high-confidence PSMs)
        /// </summary>
        private readonly HashSet<PeptideWithSetModifications> _fdrFilteredPeptides;

        private readonly List<PeptideSpectralMatch> _fdrFilteredPsms;
        private readonly List<PeptideSpectralMatch> _allPsms;
        private const double FdrCutoffForParsimony = 0.01;

        /// <summary>
        /// User-selectable option that treats differently-modified forms of a peptide as different peptides for the purposes of parsimony
        /// </summary>
        private readonly bool _treatModPeptidesAsDifferentPeptides;

        private readonly bool _useOrfInfoInProteinInference;

        public ProteinParsimonyEngine(List<PeptideSpectralMatch> allPsms, bool modPeptidesAreDifferent, CommonParameters commonParameters,
            List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters, bool useOrfInfoInProteinInference, List<string> nestedIds)
            : base(commonParameters, fileSpecificParameters, nestedIds)
        {
            _treatModPeptidesAsDifferentPeptides = modPeptidesAreDifferent;
            _useOrfInfoInProteinInference = useOrfInfoInProteinInference;

            if (!allPsms.Any())
            {
                _fdrFilteredPsms = new List<PeptideSpectralMatch>();
            }

            // parsimony will only use non-ambiguous, high-confidence PSMs
            // KEEP contaminants for use in parsimony!
            if (modPeptidesAreDifferent)
            {
                _fdrFilteredPsms = allPsms.Where(p => p.FullSequence != null && p.FdrInfo.QValue <= FdrCutoffForParsimony && p.FdrInfo.QValueNotch <= FdrCutoffForParsimony).ToList();
            }
            else
            {
                _fdrFilteredPsms = allPsms.Where(p => p.BaseSequence != null && p.FdrInfo.QValue <= FdrCutoffForParsimony && p.FdrInfo.QValueNotch <= FdrCutoffForParsimony).ToList();
            }

            // peptides to use in parsimony = peptides observed in high-confidence PSMs (including decoys)
            _fdrFilteredPeptides = new HashSet<PeptideWithSetModifications>();
            foreach (var psm in _fdrFilteredPsms)
            {
                foreach (var peptide in psm.BestMatchingPeptides.Select(p => p.Peptide))
                {
                    _fdrFilteredPeptides.Add(peptide);
                }
            }

            // we're storing all PSMs (not just FDR-filtered ones) here because we will remove some protein associations
            // from low-confidence PSMs if they can be explained by a parsimonious protein
            _allPsms = allPsms;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            ProteinParsimonyResults myAnalysisResults = new ProteinParsimonyResults(this);

            myAnalysisResults.ProteinGroups = RunProteinParsimonyEngine();

            return myAnalysisResults;
        }

        /// <summary>
        /// TODO: Summarize parsimony;
        /// Parsimony algorithm based on: https://www.ncbi.nlm.nih.gov/pubmed/14632076 Anal Chem. 2003 Sep 1;75(17):4646-58.
        /// TODO: Note describing that peptide objects with the same sequence are associated with different proteins
        /// </summary>
        private List<ProteinGroup> RunProteinParsimonyEngine()
        {
            // parsimonious list of proteins built by this protein parsimony engine
            HashSet<Protein> parsimoniousProteinList = new HashSet<Protein>();

            // if transcriptomic data is used, this is a list of proteins that would normally been removed based on parsimony algorithm,
            // but have sufficient transcriptomic support to remain
            HashSet<Protein> rescuedProteinList = new HashSet<Protein>();

            // list of peptides that can only be digestion products of one protein in the proteome (considering different protease digestion rules)
            HashSet<PeptideWithSetModifications> uniquePeptides = new HashSet<PeptideWithSetModifications>();

            // if there are no peptides observed, there are no proteins; return an empty list of protein groups
            if (_fdrFilteredPeptides.Count == 0)
            {
                return new List<ProteinGroup>();
            }

            // Parsimony stage 0: create peptide-protein associations if needed because the user wants a modification-agnostic parsimony
            // this is needed for edge cases digesting a protein .xml from UniProt that has peptide sequences shared between proteins
            // that have unevenly-shared modifications
            if (!_treatModPeptidesAsDifferentPeptides)
            {
                foreach (var protease in _fdrFilteredPsms.GroupBy(p => p.DigestionParams.Protease))
                {
                    Dictionary<string, List<PeptideSpectralMatch>> sequenceWithPsms = new Dictionary<string, List<PeptideSpectralMatch>>();

                    // for each protease, match the base sequence of each peptide to its PSMs
                    foreach (PeptideSpectralMatch psm in protease)
                    {
                        if (sequenceWithPsms.TryGetValue(psm.BaseSequence, out List<PeptideSpectralMatch> peptidesForThisBaseSequence))
                        {
                            peptidesForThisBaseSequence.Add(psm);
                        }
                        else
                        {
                            sequenceWithPsms[psm.BaseSequence] = new List<PeptideSpectralMatch> { psm };
                        }
                    }

                    var sequenceWithPsmsList = sequenceWithPsms.ToList();

                    // create new peptide-protein associations as needed
                    Parallel.ForEach(Partitioner.Create(0, sequenceWithPsmsList.Count),
                        new ParallelOptions { MaxDegreeOfParallelism = CommonParameters.MaxThreadsToUsePerFile },
                        (range, loopState) =>
                        {
                            for (int i = range.Item1; i < range.Item2; i++)
                            {
                                var baseSequence = sequenceWithPsmsList[i];

                                var peptidesWithNotchInfo = baseSequence.Value.SelectMany(p => p.BestMatchingPeptides).Distinct().ToList();

                                // if the base seq has >1 PeptideWithSetMods object and has >0 mods, it might need to be matched to new proteins
                                if (peptidesWithNotchInfo.Count > 1 && peptidesWithNotchInfo.Any(p => p.Peptide.NumMods > 0))
                                {
                                    bool needToAddPeptideToProteinAssociations = false;

                                    // numProteinsForThisBaseSequence is the total number of proteins that this base sequence is a digestion product of
                                    int numProteinsForThisBaseSequence = peptidesWithNotchInfo.Select(p => p.Peptide.Protein).Distinct().Count();

                                    if (numProteinsForThisBaseSequence == 1)
                                    {
                                        continue;
                                    }

                                    foreach (var psm in baseSequence.Value)
                                    {
                                        // numProteinsForThisPsm is the number of proteins that this PSM's peptides are associated with
                                        int numProteinsForThisPsm = psm.BestMatchingPeptides.Select(p => p.Peptide.Protein).Distinct().Count();

                                        if (numProteinsForThisPsm != numProteinsForThisBaseSequence)
                                        {
                                            // this PSM is not matched to all the proteins that it should be matched to
                                            // at this point we know that we need to make some new peptide-protein associations
                                            needToAddPeptideToProteinAssociations = true;
                                        }
                                    }

                                    if (!needToAddPeptideToProteinAssociations)
                                    {
                                        continue;
                                    }

                                    // this gets the digestion info for all of the peptide-protein associations that should exist
                                    var proteinToPeptideInfo = new Dictionary<Protein,
                                        (DigestionParams DigestParams, int OneBasedStart, int OneBasedEnd, int MissedCleavages, int Notch,
                                        CleavageSpecificity CleavageSpecificity)>();

                                    foreach (PeptideSpectralMatch psm in baseSequence.Value)
                                    {
                                        foreach (var peptideWithNotch in psm.BestMatchingPeptides)
                                        {
                                            PeptideWithSetModifications peptide = peptideWithNotch.Peptide;
                                            Protein protein = peptide.Protein;

                                            if (!proteinToPeptideInfo.ContainsKey(protein))
                                            {
                                                proteinToPeptideInfo.Add(protein,
                                                    (peptideWithNotch.Peptide.DigestionParams,
                                                    peptideWithNotch.Peptide.OneBasedStartResidueInProtein,
                                                    peptideWithNotch.Peptide.OneBasedEndResidueInProtein,
                                                    peptideWithNotch.Peptide.MissedCleavages,
                                                    peptideWithNotch.Notch,
                                                    peptideWithNotch.Peptide.CleavageSpecificityForFdrCategory));
                                            }
                                        }
                                    }

                                    // create any new associations that need to be made
                                    foreach (PeptideSpectralMatch psm in baseSequence.Value)
                                    {
                                        PeptideWithSetModifications originalPeptide = psm.BestMatchingPeptides.First().Peptide;
                                        List<MatchedFragmentIon> mfi = psm.PeptidesToMatchingFragments[originalPeptide];
                                        HashSet<Protein> psmProteins = new HashSet<Protein>(psm.BestMatchingPeptides.Select(p => p.Peptide.Protein));

                                        foreach (var proteinWithDigestInfo in proteinToPeptideInfo)
                                        {
                                            if (!psmProteins.Contains(proteinWithDigestInfo.Key))
                                            {
                                                var pep = new PeptideWithSetModifications(
                                                    proteinWithDigestInfo.Key,
                                                    proteinWithDigestInfo.Value.DigestParams,
                                                    proteinWithDigestInfo.Value.OneBasedStart,
                                                    proteinWithDigestInfo.Value.OneBasedEnd,
                                                    proteinWithDigestInfo.Value.CleavageSpecificity,
                                                    originalPeptide.PeptideDescription,
                                                    proteinWithDigestInfo.Value.MissedCleavages,
                                                    originalPeptide.AllModsOneIsNterminus,
                                                    originalPeptide.NumFixedMods);

                                                lock (_fdrFilteredPeptides)
                                                {
                                                    _fdrFilteredPeptides.Add(pep);
                                                }

                                                psm.AddProteinMatch((proteinWithDigestInfo.Value.Notch, pep), mfi);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    );
                }
            }

            // Parsimony stage 1: add proteins with unique peptides (for each protease)
            var peptidesGroupedByProtease = _fdrFilteredPeptides.GroupBy(p => p.DigestionParams.Protease);
            foreach (var peptidesForThisProtease in peptidesGroupedByProtease)
            {
                Dictionary<string, List<Protein>> peptideSequenceToProteinsForThisProtease = new Dictionary<string, List<Protein>>();
                Dictionary<string, List<PeptideWithSetModifications>> sequenceToPwsm = new Dictionary<string, List<PeptideWithSetModifications>>();

                foreach (PeptideWithSetModifications peptide in peptidesForThisProtease)
                {
                    string sequence = peptide.BaseSequence;
                    if (_treatModPeptidesAsDifferentPeptides)
                    {
                        //these and next set to full sequence but might be base sequence. treat modified as unique makes sense to use full
                        sequence = peptide.FullSequence;
                    }

                    if (peptideSequenceToProteinsForThisProtease.TryGetValue(sequence, out List<Protein> proteinsForThisPeptideSequence))
                    {
                        proteinsForThisPeptideSequence.Add(peptide.Protein);
                    }
                    else
                    {
                        peptideSequenceToProteinsForThisProtease.Add(sequence, new List<Protein> { peptide.Protein });
                    }

                    if (sequenceToPwsm.TryGetValue(sequence, out List<PeptideWithSetModifications> peptidesForThisSequence))
                    {
                        peptidesForThisSequence.Add(peptide);
                    }
                    else
                    {
                        sequenceToPwsm.Add(sequence, new List<PeptideWithSetModifications> { peptide });
                    }
                }

                foreach (var uniquePeptide in peptideSequenceToProteinsForThisProtease.Where(p => p.Value.Count == 1))
                {
                    // add the protein with the unique peptide to the parsimonious protein list
                    Protein proteinWithUniquePeptideSequence = uniquePeptide.Value.First();
                    parsimoniousProteinList.Add(proteinWithUniquePeptideSequence);

                    // add the unique peptide to the list of unique peptides
                    PeptideWithSetModifications uniquePwsm = sequenceToPwsm[uniquePeptide.Key].First();
                    uniquePeptides.Add(uniquePwsm);
                }
            }

            // Parsimony stage 2: build the peptide-protein matching structure for the parsimony greedy algorithm
            // and remove all peptides observed by proteins with unique peptides
            Dictionary<ParsimonySequence, List<Protein>> peptideSequenceToProteins = new Dictionary<ParsimonySequence, List<Protein>>();

            // this dictionary associates proteins w/ all peptide sequences (list will NOT shrink over time)
            // this is used in case of greedy algorithm ties to figure out which protein has more total peptides observed
            Dictionary<Protein, HashSet<ParsimonySequence>> proteinToPepSeqMatch = new Dictionary<Protein, HashSet<ParsimonySequence>>();

            foreach (var peptide in _fdrFilteredPeptides)
            {
                ParsimonySequence sequence = new ParsimonySequence(peptide, _treatModPeptidesAsDifferentPeptides);

                if (peptideSequenceToProteins.TryGetValue(sequence, out List<Protein> proteinsForThisPeptideSequence))
                {
                    proteinsForThisPeptideSequence.Add(peptide.Protein);
                }
                else
                {
                    peptideSequenceToProteins.Add(sequence, new List<Protein> { peptide.Protein });
                }

                if (proteinToPepSeqMatch.TryGetValue(peptide.Protein, out var peptideSequences))
                {
                    peptideSequences.Add(sequence);
                }
                else
                {
                    proteinToPepSeqMatch.Add(peptide.Protein, new HashSet<ParsimonySequence> { sequence });
                }
            }
            
            if (_useOrfInfoInProteinInference)
            {
                Dictionary<Protein, HashSet<Protein>> candidateRescueProteins = new Dictionary<Protein, HashSet<Protein>>();
                HashSet<ParsimonySequence> toRemove = new HashSet<ParsimonySequence>();
                foreach (var seq in peptideSequenceToProteins)
                {
                    bool observedAlready = seq.Value.Any(p => parsimoniousProteinList.Contains(p));

                    if (observedAlready)
                    {
                        if (seq.Value.Count > 1)
                        {
                            List<Protein> uniqueProteins = seq.Value.Where(p => parsimoniousProteinList.Contains(p)).ToList();
                            if (seq.Value.Count > uniqueProteins.Count)
                            {
                                List<Protein> sharedProteins = seq.Value.Except(uniqueProteins).ToList();
                                foreach (var uniqueProt in uniqueProteins)
                                {
                                    if (candidateRescueProteins.ContainsKey(uniqueProt))
                                    {
                                        foreach (var sharedProt in sharedProteins)
                                        {
                                            candidateRescueProteins[uniqueProt].Add(sharedProt);
                                        }
                                    }
                                    else
                                    {
                                        candidateRescueProteins.Add(uniqueProt, sharedProteins.ToHashSet());
                                    }
                                }
                            }                           

                        }
                        toRemove.Add(seq.Key);
                    }
                }
                var ProteinTranscriptAbundance = GlobalVariables.ProteinToProteogenomicInfo;
                foreach (var entry in candidateRescueProteins)
                {
                    List<Protein> toRescue = entry.Value.ToList().Where(p => ProteinTranscriptAbundance[p.Accession].CPM >= 10.0).ToList();
                    foreach (var protein in toRescue)
                    {
                        if (protein.IsDecoy == false)
                        {
                            var decoyP = proteinToPepSeqMatch.Where(p => p.Key.Accession == "DECOY_" + protein.Accession).ToList();
                            if (decoyP.Count != 0)
                            {
                                rescuedProteinList.Add(decoyP.FirstOrDefault().Key);
                            }
                            rescuedProteinList.Add(protein);
                        }
                        
                    }
                }

            }

            else
            {
                // remove the peptides observed by proteins with unique peptides
                HashSet<ParsimonySequence> toRemove = new HashSet<ParsimonySequence>();
                foreach (var seq in peptideSequenceToProteins)
                {
                    bool observedAlready = seq.Value.Any(p => parsimoniousProteinList.Contains(p));

                    if (observedAlready)
                    {
                        toRemove.Add(seq.Key);
                    }
                }
                foreach (var sequence in toRemove)
                {
                    peptideSequenceToProteins.Remove(sequence);
                }
            }           

            if (peptideSequenceToProteins.Any())
            {
                // Parsimony stage 3: greedy algorithm

                // dictionary with proteins as keys and list of associated peptide sequences as the values.
                // this data structure makes parsimony easier because the algorithm can look up a protein's peptides
                // to remove them from the list of available peptides. this list will shrink as the algorithm progresses
                var algDictionary = new Dictionary<Protein, HashSet<string>>();
                var algDictionaryProtease = new Dictionary<Protein, HashSet<ParsimonySequence>>();
                foreach (var kvp in peptideSequenceToProteins)
                {
                    foreach (var protein in kvp.Value)
                    {
                        if (algDictionaryProtease.TryGetValue(protein, out HashSet<ParsimonySequence> peptideSequencesWithProtease))
                        {
                            peptideSequencesWithProtease.Add(kvp.Key);
                        }
                        else
                        {
                            algDictionaryProtease.Add(protein, new HashSet<ParsimonySequence> { kvp.Key });
                        }

                        if (algDictionary.TryGetValue(protein, out HashSet<string> peptideSequences))
                        {
                            peptideSequences.Add(kvp.Key.Sequence);
                        }
                        else
                        {
                            algDictionary.Add(protein, new HashSet<string> { kvp.Key.Sequence });
                        }
                    }
                }

                // *** greedy algorithm loop
                int numNewSeqs = algDictionary.Max(p => p.Value.Count);
                while (numNewSeqs != 0)
                {
                    Protein bestProtein;

                    if (_useOrfInfoInProteinInference)
                    {
                        bestProtein = algDictionary
                            .OrderByDescending(p => p.Value.Count // order by num new peptide sequences multiplied by transcript weight
                            * (GlobalVariables.ProteinToProteogenomicInfo.ContainsKey(p.Key.Accession) ?
                                GlobalVariables.ProteinToProteogenomicInfo[p.Key.Accession].ProteinInferenceWeight : LongReadInfo.ProteinInferenceWeightForNoTranscriptomicsData))
                            .ThenByDescending(p => proteinToPepSeqMatch[p.Key].Count) // then order by num total sequences in protein
                            .First().Key;
                    }
                    else
                    {
                        // get the next best protein based on:
                        // 1. the number of new peptide sequences and then (in case of a tie),
                        // 2. the number of total peptides observed for the protein, regardless if they're unaccounted for or not
                        bestProtein = algDictionary
                            .Where(p => p.Value.Count == numNewSeqs)
                            .OrderByDescending(p => proteinToPepSeqMatch[p.Key].Count)
                            .First().Key;
                    }

                    parsimoniousProteinList.Add(bestProtein);

                    // remove observed peptide seqs
                    List<ParsimonySequence> observedPeptides = algDictionaryProtease[bestProtein].ToList();
                    foreach (ParsimonySequence peptide in observedPeptides)
                    {
                        List<Protein> proteinsWithThisPeptide = peptideSequenceToProteins[peptide];

                        foreach (Protein protein in proteinsWithThisPeptide)
                        {
                            algDictionary[protein].Remove(peptide.Sequence);
                            algDictionaryProtease[protein].Remove(peptide);
                        }
                    }

                    algDictionary.Remove(bestProtein);
                    algDictionaryProtease.Remove(bestProtein);

                    numNewSeqs = algDictionary.Any() ? algDictionary.Max(p => p.Value.Count) : 0;
                }

                // *** done with greedy algorithm
                parsimoniousProteinList.UnionWith(rescuedProteinList);
                // Parsimony stage 4: add back indistinguishable proteins (proteins that have identical peptide sets as parsimonious proteins)
                var allProteinsGroupedByNumPeptides = proteinToPepSeqMatch.GroupBy(p => p.Value.Count);
                var parsimonyProteinsGroupedByNumPeptides = parsimoniousProteinList.GroupBy(p => proteinToPepSeqMatch[p].Count);
                var indistinguishableProteins = new ConcurrentBag<Protein>();
                var rescueProteins = new Dictionary<Protein, ConcurrentBag<Protein>>();

                foreach (var group in allProteinsGroupedByNumPeptides)
                {
                    var parsimonyProteinsWithSameNumPeptides = parsimonyProteinsGroupedByNumPeptides.FirstOrDefault(p => p.Key == group.Key);
                    var list = group.ToList();
                    var parsimonyProteinsWithFewerPeptides = parsimonyProteinsGroupedByNumPeptides.Where(p => p.Key > group.Key).ToList();
                    if (parsimonyProteinsWithSameNumPeptides != null)
                    {
                        Parallel.ForEach(Partitioner.Create(0, list.Count),
                            new ParallelOptions { MaxDegreeOfParallelism = CommonParameters.MaxThreadsToUsePerFile },
                            (range, loopState) =>
                            {
                                for (int i = range.Item1; i < range.Item2; i++)
                                {
                                    Protein otherProtein = list[i].Key;

                                    foreach (var parsimonyProtein in parsimonyProteinsWithSameNumPeptides)
                                    {
                                        // if the two proteins have the same set of peptide sequences, they're indistinguishable
                                        // ignoring proteins that have identical gene names
                                        
                                        if (parsimonyProtein != otherProtein && proteinToPepSeqMatch[parsimonyProtein].SetEquals(proteinToPepSeqMatch[otherProtein]))
                                        {                                            
                                            indistinguishableProteins.Add(otherProtein);                                           
                                        }                                        
                                        
                                    }
                                }
                            }
                        );
                    }

                    if (parsimonyProteinsWithFewerPeptides != null)
                    {
                        Parallel.ForEach(Partitioner.Create(0, list.Count),
                            new ParallelOptions { MaxDegreeOfParallelism = CommonParameters.MaxThreadsToUsePerFile },
                            (range, loopState) =>
                            {
                                for (int i = range.Item1; i < range.Item2; i++)
                                {
                                    Protein otherProtein = list[i].Key;

                                    foreach (var peptideCount in parsimonyProteinsWithFewerPeptides)
                                    {
                                        foreach (var parsimonyProtein in peptideCount)
                                        {
                                            if (_useOrfInfoInProteinInference)
                                            {
                                                double otherProteinWeight = GlobalVariables.ProteinToProteogenomicInfo.ContainsKey(otherProtein.Accession) ?
                                                        GlobalVariables.ProteinToProteogenomicInfo[otherProtein.Accession].CPM : 0;

                                                double parsimonyProteinWeight = GlobalVariables.ProteinToProteogenomicInfo.ContainsKey(parsimonyProtein.Accession) ?
                                                    GlobalVariables.ProteinToProteogenomicInfo[parsimonyProtein.Accession].CPM : 0;

                                                var a = proteinToPepSeqMatch[otherProtein];
                                                var b = proteinToPepSeqMatch[parsimonyProtein];
                                                var sub = a.IsProperSubsetOf(b);

                                                if (proteinToPepSeqMatch[otherProtein].IsProperSubsetOf(proteinToPepSeqMatch[parsimonyProtein]))
                                                {
                                                    if (otherProteinWeight >= 10.0)
                                                    {
                                                        lock (rescueProteins)
                                                        {
                                                            if (rescueProteins.ContainsKey(parsimonyProtein))
                                                            {
                                                                rescueProteins[parsimonyProtein].Add(otherProtein);
                                                            }
                                                            else
                                                            {
                                                                rescueProteins.Add(parsimonyProtein, new ConcurrentBag<Protein>() { otherProtein });
                                                            }
                                                        }
                                                        
                                                    }
                                                }
                                            }
                                        }
                                        
                                    }
                                }
                            }
                        );
                    }
                   
                }

                foreach (Protein protein in indistinguishableProteins)
                {
                    parsimoniousProteinList.Add(protein);
                }

                foreach (var proteinCandidates in rescueProteins)
                {
                    HashSet<ParsimonySequence> coveredPeptides = new HashSet<ParsimonySequence>();
                    foreach (var protein in proteinCandidates.Value)
                    {
                        if (protein.IsDecoy == false)
                        {
                            var decoyP = proteinToPepSeqMatch.Where(p => p.Key.Accession == "DECOY_" + protein.Accession).ToList();
                            if (decoyP.Count != 0)
                            {
                                parsimoniousProteinList.Add(decoyP.FirstOrDefault().Key);
                            }

                        }
                        parsimoniousProteinList.Add(protein);
                        coveredPeptides.UnionWith(proteinToPepSeqMatch[protein]);
                    }

                    if (proteinToPepSeqMatch[proteinCandidates.Key].SetEquals(coveredPeptides))
                    {
                        double originalProteinWeight = GlobalVariables.ProteinToProteogenomicInfo.ContainsKey(proteinCandidates.Key.Accession) ?
                                                            GlobalVariables.ProteinToProteogenomicInfo[proteinCandidates.Key.Accession].CPM : 0;
                        if (originalProteinWeight < 10.0)
                        {
                            parsimoniousProteinList.Remove(proteinCandidates.Key);
                        }
                    }                    
                    
                }
            }

            // Parsimony stage 5: remove peptide objects that do not have proteins in the parsimonious list
            foreach (PeptideSpectralMatch psm in _allPsms)
            {
                // if this PSM has a protein in the parsimonious list, it removes the proteins NOT in the parsimonious list
                // otherwise, no proteins are removed (i.e., for PSMs that cannot be explained by a parsimonious protein,
                // no protein associations are removed)
                if (psm.BestMatchingPeptides.Any(p => parsimoniousProteinList.Contains(p.Peptide.Protein)))
                {
                    psm.TrimProteinMatches(parsimoniousProteinList);
                }
            }

            // construct protein groups
            List<ProteinGroup> proteinGroups = ConstructProteinGroups(uniquePeptides);

            // finished with parsimony
            return proteinGroups;
        }

        /// <summary>
        /// Builds protein groups after parsimony. Each protein group will only have 1 protein at this point.
        /// Indistinguishable protein groups will be merged during scoring for computational efficiency reasons
        /// (it's easier to tell when two protein groups are identical after they're scored)
        /// </summary>
        private List<ProteinGroup> ConstructProteinGroups(HashSet<PeptideWithSetModifications> uniquePeptides)
        {
            List<ProteinGroup> proteinGroups = new List<ProteinGroup>();
            var proteinToPeptidesMatching = new Dictionary<Protein, HashSet<PeptideWithSetModifications>>();

            foreach (var peptide in _fdrFilteredPeptides)
            {
                if (proteinToPeptidesMatching.TryGetValue(peptide.Protein, out HashSet<PeptideWithSetModifications> peptidesHere))
                {
                    peptidesHere.Add(peptide);
                }
                else
                {
                    proteinToPeptidesMatching.Add(peptide.Protein, new HashSet<PeptideWithSetModifications> { peptide });
                }
            }

            foreach (var kvp in proteinToPeptidesMatching)
            {
                var allPeptidesHere = proteinToPeptidesMatching[kvp.Key];
                var uniquePeptidesHere = new HashSet<PeptideWithSetModifications>(allPeptidesHere.Where(p => uniquePeptides.Contains(p)));

                proteinGroups.Add(new ProteinGroup(new HashSet<Protein> { kvp.Key }, allPeptidesHere, uniquePeptidesHere));
            }

            foreach (var proteinGroup in proteinGroups)
            {
                proteinGroup.AllPeptides.RemoveWhere(p => !proteinGroup.Proteins.Contains(p.Protein));
                proteinGroup.DisplayModsOnPeptides = _treatModPeptidesAsDifferentPeptides;
            }

            return proteinGroups;
        }
    }
}