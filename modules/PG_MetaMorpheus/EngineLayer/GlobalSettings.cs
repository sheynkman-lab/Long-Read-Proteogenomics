namespace EngineLayer
{
    public class GlobalSettings
    {
        public bool WriteExcelCompatibleTSVs { get; set; } = true;
        public bool UserHasAgreedToThermoRawFileReaderLicence { get; set; } = true; // TODO: change this into a CMD flag so Docker users can accept
    }
}