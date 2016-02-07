/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.info;

/**
 */
public enum InfoType {
    ProteinInfo,
    // Quality measures
    GDT_TS, DELTA_GDT_TS, GDT_HA, DELTA_GDT_HA, RMS, DELTA_RMS, RMS_HEAVY, DELTA_RMS_HEAVY, CHANGE,
    Z_GDT, Z_RMS,
    //properties
    SIZE(InfoValueType.INTEGER),
    NBL(InfoValueType.INTEGER),
    //energy
    ENERGY,
    BOND("bondEnergy"),
    ANGLE("angleEnergy"),
    PLANE("planeEnergy"),
    OO_PLANE("outOfPlaneEnergy"),
    TWO_TORSIONS,
    FLAT_RAMACH("flatRamachEnergy"),
    RAMACHANDRAN_SIDECHAIN("ramachandranSidechain"), RAMACH_STD("ramachSTD"),
    RAMACHANDRAN_CORE("ramachandranSidechainCore"),
    COOP_Z_RAMACHANDRAN("cooperativeZRamachandranSidechain"),
    COOP_Z_STD_RAMACHANDRAN("cooperativeZstdRamachandranSidechain"),
    PROPENSITY("compositePropensity"),
    COOP_Z_PROPENSITY("cooperativeZPropensity"),
    COOP_Z_STD_PROPENSITY("cooperativeZstdPropensity"),
    HYDROGEN_BONDS("hydrogenBonds"),
    HYDROGEN_BONDS_PAIRS("hydrogenBondsPairs"),
    HB_PUNISH_HOC_ANGLE("hydrogenBondsAnglesHOC"),
    HB_PUNIDH_OHN_ANGLE("hydrogenBondsAnglesOHN"),
    TETHER("tetherEnergy"), TETHER_ALL("tetherAll"),TETHER_BACKBONE("tetherBackbone"),
    INFLATE_ENERGY("inflateEnergy"), INFLATE_PER_SEGMENT("inflatePerSegment"), INFLATE_BY_SEGMENT("inflateBySegment"), INFLATE_BY_OTHER_MODEL("inflateByOtherModel"),
    BETA,
    HELIX_ATTRACTOR,
    RG, N_RG_HSS("N_RGhSS"), E_HSS_HCOIL("EhSShCoil"), E_HSS_BSS("EhSSbSS"), E_HSS_BCOIL("EhSSbCoil"), E_HSS_CSS("EhSScSS"), E_HSS_CCOIL("EhSScCoil"), HSS("hSS"), BSS("bSS"), CSS("cSS"),
    HSS_HCOIL("hSShCoil"), HSS_BSS("hSSbSS"), HSS_BCOIL("hSSbCoil"), HSS_CSS("hSScSS"), HSS_CCOIL("hSScCoil"),
    RAMACHANDRAN,
    SOLVATION("solvateEnergy"), SOLVATION_SC_POLAR("solvationSCpolar"), SOLVATION_BB_POLAR("solvationBBpolar"), SOLVATION_SC_CARBON("solvationSCcarbon"), SOLVATION_HB("solvationHB"), SOLVATION_STD("solvationSTD"),
    ATOMIC_PAIRWISE_PMF_SUMMA("atomicPairwisePMFSumma"), SUMMA_STD("summaStd"),
    COOPERATIVE_Z_SUMMA("cooperativeZSumma"), COOPERATIVE_SUMMA_POLAR("cooperativeSummaPolar"), COOPERATIVE_SUMMA_NON_POLAR("cooperativeSummaNonPolar"), COOPERATIVE_SUMMA_NEUTRAL("cooperativeSummaNeutral"), COOPERATIVE_SUMMA_POLAR_NN_OO("cooperativeSummaPolarNN_OO"), COOPERATIVE_SUMMA_POLAR_BB("cooperativeSummaPolarBb"),
    COOPERATIVE_Z_STD_SUMMA("cooperativeZstdSumma"), COOPERATIVE_STD_SUMMA_POLAR("cooperativeStdSummaPolar"), COOPERATIVE_STD_SUMMA_NON_POLAR("cooperativeStdSummaNonPolar"), COOPERATIVE_STD_SUMMA_NEUTRAL("cooperativeStdSummaNeutral"), COOPERATIVE_STD_SUMMA_POLAR_NN_OO("cooperativeStdSummaPolarNN_OO"), COOPERATIVE_STD_SUMMA_POLAR_BB("cooperativeStdSummaPolarBb"),
    TEMPLATE_ENERGY,
    EXCLUDED_VOL("excludedVolume"),
    QMOD("qMode"),
    RAPDF("samudralaEnergy"),RAPDF1("samudralaEnergy1"),
    CONSERVATION_CONTACTS8("conservationContacts8"), CONTACTS8("contacts8"),
    CONSERVATION_RG_RATIO("conservationRgRatio"),   CONSERVATION_H_RG_RATIO("conservation_H_RgRatio"),
    CONSERVATION_CONTACTS11("conservationContacts11"), CONTACTS11("contacts11"),
    CONSERVATION_CONTACTS15("conservationContacts15"), CONTACTS15("contacts15"),
    CONSERVATION_CONTACTS_HR("conservationContactsHr"), CONTACTS_HR("contactsHr"),
    CONSERVATION_RG_RATIO_HR("conservationRgRatioHr"),
    ONE("one"),
    DISTANCE_CONSTRAINTS("distanceConstraints"),
    STRETCH("stretch"),
    SECONDARY_STRUCTURE_FRACTION("secondaryStructureFraction"),
    LENGTH("length");

    public final InfoValueType valueType;
    public final String tag;

    private InfoType() {
        valueType = InfoValueType.DOUBLE;
        String s = toString();
        tag = s.toLowerCase();
    }

    private InfoType(InfoValueType valueType) {
        this.valueType = valueType;
        String s = toString();
        tag = s.charAt(0) + (s.substring(1, s.length() - 1).toLowerCase());
    }

    private InfoType(String tag) {
        valueType = InfoValueType.DOUBLE;
        this.tag = tag;
    }

    private InfoType(InfoValueType valueType, String tag) {
        this.valueType = valueType;
        this.tag = tag;
    }

    public boolean hasDoubleValue() {
        return (valueType == InfoValueType.DOUBLE);
    }

    public static boolean numeric(InfoType type) {
        return InfoValueType.numeric(type.valueType);
    }

    public static InfoType getType(String name) {
        for (InfoType infoType : values())
            if (infoType.toString().equals(name)) return infoType;
        return null;
    }
    public static InfoType getTypeByTag(String tag) {
        for (InfoType infoType : values())
            if (infoType.tag.equals(tag)) return infoType;
        throw new RuntimeException("Failed to find infoType for " + tag);
    }
}
