/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.pairwiseNonBondedTerms.CooperativeSumma;

import meshi.geometry.Distance;
import meshi.geometry.DistanceList;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomCore;
import meshi.parameters.AtomType;
import meshi.util.MeshiAttribute;

/**
 * Created by IntelliJ IDEA.
 * User: tetyanam
 * Date: 15/02/2010
 * Time: 12:16:41
 * To change this template use File | Settings | File Templates.
 */
public class CooperativeSummaProcesser {
    private static AtomType type;
    private static boolean atom1TypeIsPolarBackbone;
    private static boolean atom1TypeIsPolar ;
    private static boolean atom1TypeIsNonPolar;
    private static boolean atom1TypeIsNeutral;
    private static boolean atom1IsOxygen;
    private static boolean atom1IsNitrogen;


    private final DistanceMatrix distanceMatrix;
    //private DistanceLists nonBondedList;
//for cooperative energy terms
    protected  double[] atomEnergies;
    protected  double[] atomEnergiesPOLARS;
    protected  double[] atomEnergiesNonPOLARS;
    protected  double[] atomEnergiesNEUTRALS;
    //protected static double[] atomEnergiesPOLARSH;
    protected  double[] atomEnergiesPOLARSnoH;
    protected  double[] atomEnergiesPOLARSnoH_NNorOO;

    public double[] atomEnergies() {
        return atomEnergies;
    }

    public double[] atomEnergiesPOLARS() {
        return atomEnergiesPOLARS;
    }

    public double[] atomEnergiesNonPOLARS() {
        return atomEnergiesNonPOLARS;
    }

    public double[] atomEnergiesNEUTRALS() {
        return atomEnergiesNEUTRALS;
    }
    //public static double [] atomEnergiesPOLARSH(){return atomEnergiesPOLARSH;}

    public double[] atomEnergiesPOLARSnoH() {
        return atomEnergiesPOLARSnoH;
    }

    public double[] atomEnergiesPOLARSnoH_NNorOO() {
        return atomEnergiesPOLARSnoH_NNorOO;
    }


    public CooperativeSummaProcesser(DistanceMatrix distanceMatrix) {
        this.distanceMatrix = distanceMatrix;
        atomEnergies = new double[distanceMatrix.molecularSystem.size()];
        atomEnergiesPOLARS = new double[distanceMatrix.molecularSystem.size()];
        //atomEnergiesPOLARSH = new double[distanceMatrix.molecularSystem.size()];
        atomEnergiesPOLARSnoH = new double[distanceMatrix.molecularSystem.size()];
        atomEnergiesPOLARSnoH_NNorOO = new double[distanceMatrix.molecularSystem.size()];

        atomEnergiesNonPOLARS = new double[distanceMatrix.molecularSystem.size()];
        atomEnergiesNEUTRALS = new double[distanceMatrix.molecularSystem.size()];

    }

    public void reset() {
        for (int i = 0; i < atomEnergiesPOLARS.length; i++) {
            atomEnergies[i] = 0;
            atomEnergiesPOLARS[i] = 0;
            //       atomEnergiesPOLARSH[i] = 0;
            atomEnergiesPOLARSnoH[i] = 0;
            atomEnergiesPOLARSnoH_NNorOO[i] = 0;

            atomEnergiesNonPOLARS[i] = 0;
            atomEnergiesNEUTRALS[i] = 0;
        }
        /*
DistanceLists nonBondedList = distanceMatrix.nonBondedList();
for (Distance distance:nonBondedList) {
if (distance.mode().frozen) continue;
summaAttribute = (SummaAttribute) distance.getAttribute(MeshiAttribute.SUMMA_ATTRIBUTE);
if (summaAttribute != null)
                 summaAttribute.restart();
}
        */
    }

    public void setCooperativeStatistic() {

/*
//For statistic

        atomEnergies[atom1.number] += halfEnergy0;
        atomEnergies[atom2.number] += halfEnergy0;
     //for atom1
     //  if (atom2.type().isPolar() || atom2.type().isPolarSideChains()) {
          if (atom2.type().isPolar() ) {
           atomEnergiesPOLARS[atom1.number] += halfEnergy0;
       }
       else
         if (atom2.type().isNonPolar()) atomEnergiesNonPOLARS[atom1.number] += halfEnergy0;
         else
             if (atom2.type().isNeutral()) atomEnergiesNEUTRALS[atom1.number] += halfEnergy0;

     //for atom2
       if (atom1.type().isPolar() || atom1.type().isPolarSideChains()) {
         atomEnergiesPOLARS[atom2.number] += halfEnergy0;
     }
     else
       if (atom1.type().isNonPolar()) atomEnergiesNonPOLARS[atom2.number] += halfEnergy0;
       else
           if (atom1.type().isNeutral()) atomEnergiesNEUTRALS[atom2.number] += halfEnergy0;
          //*/

//--------------------------------------------------For statistic with H bond--------------------------------------------------
        SummaAttribute summaAttribute;
        int atom1Number;
        AtomType atom1Type;
        int atom2Number;
        AtomType atom2Type;
        int atom1ResidueNumber;
        int atom2ResidueNumber;
        double halfEnergy0;
        double atom1AtomEnergies,atom1AtomEnergiesPOLARS,atom1AtomEnergiesPOLARSH,atom1AtomEnergiesPOLARSnoH,atom1AtomEnergiesPOLARSnoH_NNorOO,atom1AtomEnergiesNonPOLARS,atom1AtomEnergiesNEUTRALS;

        reset();
        for (DistanceList row : distanceMatrix.nonBondedList()){
            AtomCore atom1 = row.atomOne;
            atom1Type = atom1.type();
            if (atom1Type.isHydrogen()) continue;
            atom1AtomEnergies = atom1AtomEnergiesPOLARS = atom1AtomEnergiesPOLARSH = atom1AtomEnergiesPOLARSnoH = 0;
            atom1AtomEnergiesPOLARSnoH_NNorOO = atom1AtomEnergiesNonPOLARS = atom1AtomEnergiesNEUTRALS = 0;
            atom1Number = atom1.number;
            atom1ResidueNumber = atom1.atom.residueNumber();

            atom1TypeIsPolarBackbone = atom1Type.isPolarBackbone();
            atom1TypeIsPolar         = atom1Type.isPolar();
            atom1TypeIsNonPolar      = atom1Type.isNonPolar();
            atom1TypeIsNeutral       = atom1Type.isNeutral();
            atom1IsOxygen            = atom1Type.isOxygen();
            atom1IsNitrogen          = atom1Type.isNitrogen();

            for (Distance distance : row) {
                summaAttribute = (SummaAttribute) distance.getAttribute(MeshiAttribute.SUMMA_ATTRIBUTE);
                if (summaAttribute == null) continue;
                Atom atom2 = distance.atom2();
                atom2Type =   atom2.type();
                atom2Number = atom2.number();
                halfEnergy0 = summaAttribute.halfEnergy;
                 atom2ResidueNumber = atom2.residueNumber();
                boolean atom2TypeIsPolarBackbone = atom2Type.isPolarBackbone();
                boolean atom2TypeIsPolar = atom2Type.isPolar();
                boolean atom2TypeIsNonPolar = atom2Type.isNonPolar();
                boolean atom2TypeIsNeutral = atom2Type.isNeutral();
                boolean atom2IsOxygen = atom2Type.isOxygen();
                boolean atom2IsNitrogen = atom2Type.isNitrogen();
                atom1AtomEnergies += halfEnergy0;
                atomEnergies[atom2Number] += halfEnergy0;

                if ((!atom1TypeIsPolarBackbone) || (!atom2TypeIsPolarBackbone)) {
                    if (atom2TypeIsPolar) {
                        atom1AtomEnergiesPOLARS += halfEnergy0;
                    }
                    else if (atom2TypeIsNonPolar)
                        atom1AtomEnergiesNonPOLARS += halfEnergy0;
                    else if (atom2TypeIsNeutral)
                        atom1AtomEnergiesNEUTRALS += halfEnergy0;
                    else
                        throw new RuntimeException("The weird case of polarity assignment (see AtomType.java) in " + this);
                    //for atom2
                    if (atom1TypeIsPolar)
                        atomEnergiesPOLARS[atom2Number] += halfEnergy0;
                    else if (atom1TypeIsNonPolar)
                        atomEnergiesNonPOLARS[atom2Number] += halfEnergy0;
                    else if (atom1TypeIsNeutral)
                        atomEnergiesNEUTRALS[atom2Number] += halfEnergy0;
                    else
                        throw new RuntimeException("The weird case of polarity assignment (see AtomType.java) in " + this + " #1 \n" + distance.atom1 + " " + distance.atom1.type() + "\n" + distance.atom2 + " " + distance.atom2.type());
                } else {
                    if (!(atom1TypeIsPolarBackbone && atom2TypeIsPolarBackbone))
                        throw new RuntimeException("The weird case of polarity assignment (see AtomType.java) in " + this + " #1 \n" + distance.atom1 + " " + distance.atom1.type() + "\n" + distance.atom2 + " " + distance.atom2.type());
//only two Polars N or O are stayed

                    if ((atom1IsOxygen && atom2IsOxygen)
                            ||
                            (atom1IsNitrogen && atom2IsNitrogen)) {
                        atom1AtomEnergiesPOLARSnoH_NNorOO += halfEnergy0;
                        atomEnergiesPOLARSnoH_NNorOO[atom2Number] += halfEnergy0;
                    } else {  //hBond
//                  atomEnergiesPOLARS[atom1.number] += halfEnergy0;
                        //              atomEnergiesPOLARS[atom2.number] += halfEnergy0;

                        if (Math.abs(atom1ResidueNumber - atom2ResidueNumber) < 3)
                            continue;
                        atom1AtomEnergiesPOLARSnoH += halfEnergy0;
                        atomEnergiesPOLARSnoH[atom2Number] += halfEnergy0;

                        if (!((atom1IsNitrogen && atom2IsOxygen) ||
                                (atom2IsNitrogen && atom1IsOxygen)))
                            throw new RuntimeException("The weird case of polarity assignment (see AtomType.java) in " + this + "\n" +
                                    distance + "\n" +
                                    atom1Type + "  " + atom2Type);
                    }
                }
            }
            atomEnergies[atom1Number] += atom1AtomEnergies;
            atomEnergiesPOLARS[atom1Number] = atom1AtomEnergiesPOLARS;
            atomEnergiesPOLARSnoH[atom1Number] = atom1AtomEnergiesPOLARSnoH;
            atomEnergiesPOLARSnoH_NNorOO[atom1Number] = atom1AtomEnergiesPOLARSnoH_NNorOO;
            atomEnergiesNonPOLARS[atom1Number] = atom1AtomEnergiesNonPOLARS;
            atomEnergiesNEUTRALS[atom1Number] = atom1AtomEnergiesNEUTRALS;
        }
    }

}
