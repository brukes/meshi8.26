/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.parameters;

import meshi.util.filters.Filter;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomCore;

public enum AtomType {
    //Termini 0-2
    TRN(BBatom.N, Element.N, Polar.OTHER), // Charged Nitrogen at tne N-terminus
    TRC(BBatom.C, Polar.OTHER), // Carboxyl Carbon at the C-terminus
    TRO(BBatom.O, Polar.OTHER), // Charged Oxygen at tne C-terminus
    //Ala 3-8
    AH(BBatom.H, Element.H), AN(BBatom.N, Element.N), ACA(BBatom.CA), AC(BBatom.C), AO(BBatom.O, Element.O), ACB(BBatom.CB),
    //Cys 9-15
    CH(BBatom.H, Element.H, Polar.OTHER), CN(BBatom.N, Element.N, Polar.OTHER), CCA(BBatom.CA, Polar.OTHER), CC(BBatom.C, Polar.OTHER), CO(BBatom.O, Element.O, Polar.OTHER), CCB(BBatom.CB, Polar.OTHER),
    CSG(Element.S, Polar.OTHER), // True for non-Disulfides
    //Asp 16-23
    DH(BBatom.H, Element.H), DN(BBatom.N, Element.N), DCA(BBatom.CA, Polar.NEUTRAL), DC(BBatom.C, Polar.NEUTRAL), DO(BBatom.O, Element.O), DCB(BBatom.CB, Polar.NEUTRAL),
    DCG(Polar.NEUTRAL),
    DOD(Element.O, Polar.POLAR_SIDECHAINS, Charged.NEGATIVE),
    //Glu 24-32
    EH(BBatom.H, Element.H), EN(BBatom.N, Element.N), ECA(BBatom.CA, Polar.NEUTRAL), EC(BBatom.C, Polar.NEUTRAL), EO(BBatom.O, Element.O), ECB(BBatom.CB, Polar.NEUTRAL),
    ECG(Polar.NEUTRAL),
    ECD(Polar.NEUTRAL),
    EOE(Element.O, Polar.POLAR_SIDECHAINS, Charged.NEGATIVE),
    //Phe 33-42
    FH(BBatom.H, Element.H), FN(BBatom.N, Element.N), FCA(BBatom.CA), FC(BBatom.C), FO(BBatom.O, Element.O), FCB(BBatom.CB),
    FCG,
    FCD,
    FCE,
    FCZ,
    //Gly 43-47
    GH(BBatom.H, Element.H), GN(BBatom.N, Element.N), GCA(BBatom.CA, Polar.NEUTRAL), GC(BBatom.C, Polar.NEUTRAL), GO(BBatom.O, Element.O),
    //His 48-60
    HH(BBatom.H, Element.H), HN(BBatom.N, Element.N), HCA(BBatom.CA, Polar.NEUTRAL), HC(BBatom.C, Polar.NEUTRAL), HO(BBatom.O, Element.O), HCB(BBatom.CB, Polar.NEUTRAL),
    HCG(Polar.NEUTRAL),
    HCD(Polar.NEUTRAL),
    HND(Element.N, Polar.POLAR_SIDECHAINS),
    HHD(Element.H),
    HCE(Polar.NEUTRAL),
    HNE(Element.N, Polar.POLAR_SIDECHAINS),
    HHE(Element.H),
    //Ile 61-69
    IH(BBatom.H, Element.H), IN(BBatom.N, Element.N), ICA(BBatom.CA), IC(BBatom.C), IO(BBatom.O, Element.O), ICB(BBatom.CB),
    ICG1,
    ICG2,
    ICD,
    //Lys 70-79
    KH(BBatom.H, Element.H), KN(BBatom.N, Element.N), KCA(BBatom.CA, Polar.NEUTRAL), KC(BBatom.C, Polar.NEUTRAL), KO(BBatom.O, Element.O), KCB(BBatom.CB, Polar.NEUTRAL),
    KCG(Polar.NEUTRAL),
    KCD(Polar.NEUTRAL),
    KCE(Polar.NEUTRAL),
    KNZ(Element.N, Polar.POLAR_SIDECHAINS, Charged.POSITIVE),
    //Leu 80-88
    LH(BBatom.H, Element.H), LN(BBatom.N, Element.N), LCA(BBatom.CA), LC(BBatom.C), LO(BBatom.O, Element.O), LCB(BBatom.CB),
    LCG,
    LCD1,
    LCD2,
    //Met 89-97
    MH(BBatom.H, Element.H), MN(BBatom.N, Element.N), MCA(BBatom.CA), MC(BBatom.C), MO(BBatom.O, Element.O), MCB(BBatom.CB),
    MCG,
    MSD(Element.S, Polar.NONPOLAR),
    MCE,
    //Asn 98-108
    NH(BBatom.H, Element.H), NN(BBatom.N, Element.N), NCA(BBatom.CA, Polar.NEUTRAL), NC(BBatom.C, Polar.NEUTRAL), NO(BBatom.O, Element.O), NCB(BBatom.CB, Polar.NEUTRAL),
    NCG(Polar.NEUTRAL),
    NOD(Element.O, Polar.POLAR_SIDECHAINS),
    NND(Element.N, Polar.POLAR_SIDECHAINS),
    NHD1(Element.H),
    NHD2(Element.H),
    //Pro 109-115
    PN(BBatom.N, Element.N, Polar.POLAR_SIDECHAINS), PCA(BBatom.CA, Polar.NEUTRAL), PC(BBatom.C, Polar.NEUTRAL), PO(BBatom.O, Element.O, Polar.POLAR_SIDECHAINS), PCB(BBatom.CB, Polar.NEUTRAL),
    PCG(Polar.NEUTRAL),
    PCD(Polar.NEUTRAL),

    //Gln 116-127
    QH(BBatom.H, Element.H), QN(BBatom.N, Element.N), QCA(BBatom.CA, Polar.NEUTRAL), QC(BBatom.C, Polar.NEUTRAL), QO(BBatom.O, Element.O), QCB(BBatom.CB, Polar.NEUTRAL),
    QCG(Polar.NEUTRAL),
    QCD(Polar.NEUTRAL),
    QOE(Element.O, Polar.POLAR_SIDECHAINS),
    QNE(Element.N, Polar.POLAR_SIDECHAINS),
    QHE1(Element.H),
    QHE2(Element.H),

    //Arg 128-139
    RH(BBatom.H, Element.H), RN(BBatom.N, Element.N), RCA(BBatom.CA, Polar.NEUTRAL), RC(BBatom.C, Polar.NEUTRAL), RO(BBatom.O, Element.O), RCB(BBatom.CB, Polar.NEUTRAL),
    RCG(Polar.NEUTRAL),
    RCD(Polar.NEUTRAL),
    RNE(Element.N, Polar.POLAR_SIDECHAINS, Charged.POSITIVE),
    RHE(Element.H),
    RCZ(Polar.NEUTRAL),
    RNH(Element.N, Polar.POLAR_SIDECHAINS, Charged.POSITIVE),
    //Ser
    SH(BBatom.H, Element.H), SN(BBatom.N, Element.N), SCA(BBatom.CA, Polar.NEUTRAL), SC(BBatom.C, Polar.NEUTRAL), SO(BBatom.O, Element.O), SCB(BBatom.CB, Polar.NEUTRAL),
    SOG(Element.O, Polar.POLAR_SIDECHAINS),
    //Thr
    TH(BBatom.H, Element.H), TN(BBatom.N, Element.N), TCA(BBatom.CA, Polar.NEUTRAL), TC(BBatom.C, Polar.NEUTRAL), TO(BBatom.O, Element.O), TCB(BBatom.CB, Polar.NEUTRAL),
    TCG(Polar.NEUTRAL),
    TOG(Element.O, Polar.POLAR_SIDECHAINS),
    //Val
    VH(BBatom.H, Element.H), VN(BBatom.N, Element.N), VCA(BBatom.CA), VC(BBatom.C), VO(BBatom.O, Element.O), VCB(BBatom.CB),
    VCG1,
    VCG2,
    //Trp
    WH(BBatom.H, Element.H), WN(BBatom.N, Element.N), WCA(BBatom.CA), WC(BBatom.C), WO(BBatom.O, Element.O), WCB(BBatom.CB),
    WCG,
    WCD1(),
    WCD2,
    WCE2(),
    WCE3,
    WNE(Element.N, Polar.POLAR_SIDECHAINS),//Polar.POLAR),
    WHE(Element.H),
    WCZ2,
    WCZ3,
    WCH2,
    //Tyr
    YH(BBatom.H, Element.H), YN(BBatom.N, Element.N), YCA(BBatom.CA), YC(BBatom.C), YO(BBatom.O, Element.O), YCB(BBatom.CB),
    YCG,
    YCD,
    YCE,
    YCZ(),
    YOH(Element.O, Polar.POLAR_SIDECHAINS),
    // generic
    XXX,
    // for helium cluster minimization
    HE;

    private final BBatom bbAtom;
    private final Element element;
    private final Polar polar;
    private final Charged charged;

    //------------------- new constructors -------------------

    // constructors
//full constructor
    private AtomType(BBatom bbAtom, Element element, Polar polar, Charged charged) {
        this.bbAtom = bbAtom;
        this.element = element;
        this.polar = polar;
        this.charged = charged;
    }

    private AtomType(Element element, Polar polar, Charged charged) {
        this(BBatom.NOT_BACKBONE, element, polar, charged);
    }

    // without Charged
    private AtomType(BBatom bbAtom, Element element, Polar polar) {
        this(bbAtom, element, polar, Charged.NEUTRAL);
    }

    private AtomType(BBatom bbAtom, Element element) {
        this(bbAtom, element, getPolar(element));
    }

    private AtomType(Element element) {
        this(BBatom.NOT_BACKBONE, element);
    }

    private AtomType(BBatom bbAtom) {
        this(bbAtom, Element.C);
    }

    private AtomType(Element element, Charged charged) {
        this(BBatom.NOT_BACKBONE, element, getPolar(element), charged);
    }

    private AtomType(Polar polar) {
        this(BBatom.NOT_BACKBONE, Element.C, polar);
    }

    private AtomType(BBatom bbAtom, Polar polar) {
        this(bbAtom, Element.C, polar);
    }

    private AtomType() {
        this(BBatom.NOT_BACKBONE);
    }

    private AtomType(Element element, Polar polar) {
        this(BBatom.NOT_BACKBONE, element, polar);
    }

    // A helping method for the constructors
    private static Polar getPolar(Element element) {
        Polar polar;
        switch (element) {
            case H:
                polar = Polar.OTHER;
                break;
            case N:
            case O:
                polar = Polar.POLAR;
                break;
            default:
                polar = Polar.NONPOLAR;
        }
        return polar;
    }
    //-------------------------------------------------------------


    //
    public final boolean backbone() {
        return bbAtom != BBatom.NOT_BACKBONE;
    }

    public final boolean backboneH() {
        return bbAtom == BBatom.H;
    }

    public final boolean backboneN() {
        return bbAtom == BBatom.N;
    }

    public final boolean backboneCA() {
        return bbAtom == BBatom.CA;
    }

    public final boolean backboneC() {
        return bbAtom == BBatom.C;
    }

    public final boolean backboneO() {
        return bbAtom == BBatom.O;
    }

    public final boolean isCarbon() {
        return element == Element.C;
    }

    public final boolean isOxygen() {
        return element == Element.O;
    }

    public final boolean isNitrogen() {
        return element == Element.N;
    }

    public final boolean isSulfur() {
        return element == Element.S;
    }

    public final boolean isHydrogen() {
        return element == Element.H;
    }

    public final BBatom bbAtom() {
        return bbAtom;
    }

    public final boolean isPolar() {
        return ((polar == Polar.POLAR) || (polar == Polar.POLAR_SIDECHAINS));
    }

    public final boolean isPolarBackbone() {
        return polar == Polar.POLAR;
    }

    public final boolean isPolarSideChains() {
        return polar == Polar.POLAR_SIDECHAINS;
    }

    public final boolean isNonPolar() {
        return polar == Polar.NONPOLAR;
    }

    public final boolean isNeutral() {
        return polar == Polar.NEUTRAL;
    }

    public final boolean isOther() {
        return polar == Polar.OTHER;
    }

    public final boolean isCharged() {
        return charged == Charged.POSITIVE || charged == Charged.NEGATIVE;
    }


    public static int numberOfHeavyAtomTypes() {
        int counter = 0;
        for (AtomType type : values())
            if (! type.isHydrogen()) counter++;
        return counter;
    }
    public static AtomType type(String residueName, String atomName) {
        if (atomName.equals("OXT")) return TRO;
        ResidueType residueType = ResidueType.type(residueName);
        String name = residueType.nameOneLetter() + atomName;
        for (AtomType type : AtomType.values()) {
            if (name.equals(type.toString())) return type;
            else if ((name.length() >= 4) &&
                    name.substring(0, 3).equals(type.toString())) return type;
        }
        if (residueType == ResidueType.ASN) {
            if (atomName.equals("1HD2")) return AtomType.NHD1;
            if (atomName.equals("HD21")) return AtomType.NHD1;
            if (atomName.equals("2HD2")) return AtomType.NHD2;
            if (atomName.equals("HD22")) return AtomType.NHD2;
        }
        if (residueType == ResidueType.GLN) {
            if (atomName.equals("1HE2")) return AtomType.QHE1;
            if (atomName.equals("HE21")) return AtomType.QHE1;
            if (atomName.equals("2HE2")) return AtomType.QHE2;
            if (atomName.equals("HE22")) return AtomType.QHE2;
        }
        return XXX;
    }

    public static AtomType type(String name) {
        for (AtomType type : AtomType.values()) {
            if (name.equals(type.toString())) return type;
        }
        return XXX;
    }

    public static AtomType type(int ordinal) {
        for (AtomType type : AtomType.values())
            if (ordinal == type.ordinal()) return type;
        return XXX;
    }

    public static class ChargedFilter implements Filter {
        public boolean accept(Object o) {
            if (o instanceof Atom) return (((Atom) o).type().isCharged());
            if (o instanceof AtomCore) return (((AtomCore) o).type().isCharged());
            throw new RuntimeException("Weird parameter to ChargedFilter.accept " + o);
        }
    }

}
//
