﻿using Proteomics;
using System;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer
{
    public class LongReadInfo
    {
        public readonly Protein Protein;
        public readonly int NumReads;

        public LongReadInfo(Protein protein, int numReads)
        {
            this.Protein = protein;
            this.NumReads = numReads;
        }
    }
}