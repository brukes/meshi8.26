
package meshi.energy.hydrogenBondPairPatterns;

import meshi.energy.hydrogenBondsPairs.HydrogenBondsPairsCreator;
import meshi.energy.hydrogenBondsPairs.HydrogenBondsPairsEnergy;
import meshi.energy.hydrogenBondsPairs.PairOfHydrogenBondsElements;
import meshi.geometry.Distance;

import java.io.*;

/**
 * Created by ofer on 31/10/15.
 */
public class HydrogenBondPairTypes {
    //-----------------data--------------------
    private String [] proteinList;
    private String [] listOfHBPTypes;
    private String [][] typeLocationsInProteins;
    private final int largestProteinSize=999;


    //------------------constructors-----------

    public HydrogenBondPairTypes (String [] proteinList) {
        this.proteinList = proteinList;
        listOfHBPTypes = listOfHBPTypesMaker();
        typeLocationsInProteins = new String [proteinList.length][listOfHBPTypes.length];

    }


    //-----------------methods-----------------

    private static String[] listOfHBPTypesMaker () {

        String[] listOfHBPTypes = new String[750];
        try {
            BufferedReader br = new BufferedReader(new FileReader("HBBetaHelixParametersSmall.txt"));
            String line = br.readLine();
            int lineNum = 0;

            while (line != null) {
                listOfHBPTypes[lineNum] = HBPTypeNamer(line, lineNum);
                line = br.readLine();

                lineNum++;
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
        return listOfHBPTypes;
    }

    private static String HBPTypeNamer (String line, int lineNum) {
        String name;
        String firstBond= ".FB".concat(line.substring(12, 15).replaceAll(" ", ""));
        line=removeSpace(line.substring(16));
        String secondBond = ".SB".concat(line.substring(13, 16).replaceAll(" ", ""));
        line=removeSpace(line.substring(17));
        String hDis= ".HD".concat(line.substring(8, 11).replaceAll(" ", ""));
        line=removeSpace(line.substring(12));
        String oDis= ".OD".concat(line.substring(8, 11).replaceAll(" ", ""));
        line=removeSpace(line.substring(12));
        String hoDis= ".HO".concat(line.substring(9, 12).replaceAll(" ", ""));
        line=removeSpace(line.substring(13));
        String ohDis= ".OH".concat(line.substring(9, 11).replaceAll(" ", ""));
        String value =".V".concat(line.substring(12).replaceAll(" ", "").replaceAll("[*]", ""));
        name=firstBond.concat(secondBond.concat(hDis.concat(oDis.concat(hoDis.concat(ohDis.concat(value))))));

        if (lineNum<671){
            name="BETA".concat(name);
        }
        else name="HELIX".concat(name);

        return name;
    }

    private static String HBPTypeNamer (PairOfHydrogenBondsElements element){
        String name;
        String secondaryStructure;

        boolean isBETA= element.isBetaPair();
        if (isBETA){
            secondaryStructure = "BETA";
        }
        else{
            secondaryStructure = "HELIX";
        }
        String test= ""+element.seqDistance1();
        name = secondaryStructure +".FB"+element.seqDistance1()+".SB"+element.seqDistance2()+".HD"+element.hhSeqDistance()+
                ".OD"+element.ooSeqDistance()+".HO"+element.hoSeqDistance()+".OH"+element.ohSeqDistance()+".V"+element.getPairValue();

        return name;
    }

    private static String removeSpace (String line) {
        while (line.charAt(0) == ' ') {
            line = line.substring(1);
        }
        return line;
    }

    public static void insertTypeLocation (HydrogenBondsPairsCreator hydrogenBondsPairsCreator, String pdbName){

    }

    private void patternStringMatrixUpdater(HydrogenBondsPairsCreator hydrogenBondsPairsCreator, String pdbName)
            throws FileNotFoundException, UnsupportedEncodingException {

        int i = 0;
        int matrixI = findProteinIndex(pdbName);
        for (PairOfHydrogenBondsElements elementT : ((HydrogenBondsPairsEnergy) hydrogenBondsPairsCreator.term()).getPairsOfHBEElementsList()) {
            Distance hbt1 = elementT.HOelement1();
            Distance hbt2 = elementT.HOelement2();
            if (elementT.isGoodPair()) {  // added this if statement isGoodPair
                int j = 0;
                for (PairOfHydrogenBondsElements elementS : ((HydrogenBondsPairsEnergy) hydrogenBondsPairsCreator.term()).getPairsOfHBEElementsList()) {
                    Distance hbs1 = elementS.HOelement1();
                    Distance hbs2 = elementS.HOelement2();
                    if (i < j && elementS.isGoodPair() && (hbt1.equals(hbs1) || hbt1.equals(hbs2) || hbt2.equals(hbs1) || hbt2.equals(hbs2))) {  //i<j was i!=j added isGoodPair
                        //find is beta boolean
                        //String elementIName=HBPTypeNamer(elementI);

                        int matrixJ1 = findPatternIndex(elementT);
                        int matrixJ2 = findPatternIndex(elementS);
                        if (matrixI != -1 && matrixJ1 != -1 && matrixJ2 !=-1) {
                            if (typeLocationsInProteins [matrixI][matrixJ1]==null){
                                typeLocationsInProteins [matrixI][matrixJ1]= hbt1.atom1Number+","+ hbt1.atom2Number()+","+
                                        hbt2.atom1Number+","+hbt2.atom2Number();
                            }
                            else {
                                typeLocationsInProteins [matrixI][matrixJ1]=typeLocationsInProteins [matrixI][matrixJ1]+","+
                                        hbt1.atom1Number+","+ hbt1.atom2Number()+","+
                                        hbt2.atom1Number+","+hbt2.atom2Number();
                            }
                            if (typeLocationsInProteins [matrixI][matrixJ2]==null){
                                typeLocationsInProteins [matrixI][matrixJ2]=hbs1.atom1Number+","+hbs1.atom2Number()+","+
                                        hbs2.atom1Number+","+hbs2.atom2Number();
                            }
                            else {
                                typeLocationsInProteins [matrixI][matrixJ2]=typeLocationsInProteins [matrixI][matrixJ2]+","+
                                        hbs1.atom1Number+","+hbs1.atom2Number()+","+
                                        hbs2.atom1Number+","+hbs2.atom2Number();
                            }
                        }
                    }
                    j++;
                }
            }
            i++;
        }
        sortLocations(matrixI);
    }

    private static int findProteinIndex (String pdbName){
        int i=-1;

        return i;
    }

    private static int findPatternIndex (PairOfHydrogenBondsElements elementI){
        int i=-1;

        return i;
    }

    private void sortLocations (int proteinIndex){
        //This function will go over all indexes for the given protein and sort the numbered locations for each HBP and remove duplicate locations

    }

}
