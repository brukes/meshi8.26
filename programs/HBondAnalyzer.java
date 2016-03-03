package programs;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.EvaluationException;
import meshi.energy.TotalEnergy;
import meshi.energy.hydrogenBond.HB_DistanceAttribute;
import meshi.energy.hydrogenBond.HydrogenBondsCreator;
import meshi.energy.hydrogenBond.HydrogenBondsEnergy;
import meshi.energy.hydrogenBondsAngle.HBondsPunishOHNAngleCreator;
import meshi.energy.hydrogenBondsAngle.HbondsPunishHOCAngleCreator;
import meshi.energy.hydrogenBondsPairs.HydrogenBondsPairsCreator;
import meshi.energy.hydrogenBondsPairs.HydrogenBondsPairsEnergy;
import meshi.energy.hydrogenBondsPairs.PairOfHydrogenBondsElements;
import meshi.energy.hydrogenBondsPairs.PairsOfHBEElementsList;
import meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma.AtomicPairwisePMFSummaCreator;
import meshi.geometry.Distance;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtoms;
import meshi.sequences.AlignmentException;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.util.*;

import java.io.*;


/**
 * Created by chen on 11/12/2014.
 */
public class HBondAnalyzer {
    public final static int NUMBRER_OF_PAIR_TYPES = 750;
    private static Protein nativeStructure = null;
    public static void main(String[] args) throws UpdateableException, EvaluationException, AlignmentException, FileNotFoundException, UnsupportedEncodingException, RuntimeException {
        MeshiProgram.initRandom(0);
        long time    = System.currentTimeMillis();
        File dsspDir = new File("nativeStructures/dssp");
        File pdbDir  = new File("nativeStructures/pdb");
        File[] dsspDirectoryListing = dsspDir.listFiles();
        File[] pdbDirectoryListing  = pdbDir.listFiles();

        System.out.println(dsspDir);
        //These matrices represent the Graph with each row/column representing a node and a non zero cell represent an edge.
        int[][] patternMatrix = new int[NUMBRER_OF_PAIR_TYPES][NUMBRER_OF_PAIR_TYPES];
        String[][] stringMatrix = new String[NUMBRER_OF_PAIR_TYPES][NUMBRER_OF_PAIR_TYPES];
        String[] pairsInProteins = new String[pdbDirectoryListing.length];
        String[] pairsLine = {""};

        for (int i=0; i<dsspDirectoryListing.length; i++){
            System.out.println("pdb "+pdbDirectoryListing[i].getName());
            System.out.println("dssp "+dsspDirectoryListing[i].getName());
        }


        if (dsspDirectoryListing != null && pdbDirectoryListing != null && dsspDirectoryListing.length==pdbDirectoryListing.length) {
            int numberOfProteins = dsspDirectoryListing.length;
            for (int iProtein=0; iProtein< numberOfProteins; iProtein++) {
                String pdbFileName = pdbDirectoryListing[iProtein].getName();
                if (!pdbFileName.startsWith(".")) {
                    String dsspFileName = dsspDir + "/" + pdbFileName + ".dssp";
                    pdbFileName = pdbDir + "/" + pdbFileName;
                    String[] arguments = {pdbFileName, dsspFileName, "commands.MCM"};

                    System.out.println("Analyzing "+iProtein+" "+pdbFileName) ;
                    run(arguments, patternMatrix, stringMatrix, pairsLine);
                    pairsInProteins[iProtein]=pairsLine[0];
                    pairsLine[0]="";
                    System.out.println("Done with "+iProtein+" "+pdbFileName) ;
                }
            }
        }
        else if (dsspDirectoryListing == null || pdbDirectoryListing == null) {
            throw new FileNotFoundException("Problem with protein folders : null listFiles");
        }
        else {
            throw new FileNotFoundException("Problem with protein folders : listFiles of different length");
        }


        String[] patternNames = patternNamesArrayMaker();
        graphMaker(patternNames, patternMatrix, stringMatrix);
        patternMatrixToFile (patternMatrix);
        graphMakerNoLocations(patternNames, patternMatrix, stringMatrix);
        locationGraphMaker(patternMatrix, stringMatrix);
        pairsInProteinsMaker(pairsInProteins);
        System.out.println("Done!");
        System.out.println("This took: "+((System.currentTimeMillis()-time)/1000)+" seconds");
    }

    public static void run (String[] args, int[][] patternMatrix, String[][] stringMatrix, String[] pairsLine) throws UpdateableException, EvaluationException, AlignmentException, FileNotFoundException, UnsupportedEncodingException{
        Protein model = new Protein(new AtomList(args[0]), ResidueExtendedAtoms.creator);
        if (nativeStructure == null) {
            String nativeStructureName = model.name().substring(0, 8) + ".N.pdb";
            File nativeStructureFile = new File(nativeStructureName);
            if (nativeStructureFile.exists()) {
                nativeStructure = new Protein(new AtomList(nativeStructureName), ResidueExtendedAtoms.creator);
            }
        }
        if (nativeStructure != null) {
            Utils.printDebug("run ", "xxxxxxxxxxxxx " + nativeStructure);
            ModelAnalyzer analyzer = new ModelAnalyzer(model, nativeStructure, model, null, ResidueAlignmentMethod.IDENTITY);
            double[] gdt = Rms.gdt(nativeStructure,model);
            Utils.printDebug("run ", "xxxxxxxxxxxxx " + nativeStructure + " " + gdt[0]);
        }

        System.out.println("new model " + model);
        Utils.AssignDSSP(model, args[1]);
        model.chain().print();
        CommandList commands = new CommandList(args[2]);
        HydrogenBondsCreator hydrogenBondsCreator = new HydrogenBondsCreator();
        HydrogenBondsPairsCreator hydrogenBondsPairsCreator = new HydrogenBondsPairsCreator(hydrogenBondsCreator);
        AtomicPairwisePMFSummaCreator summaCreator = new AtomicPairwisePMFSummaCreator();

        EnergyCreator[] energyCreators = {
                summaCreator,
                hydrogenBondsCreator,
                hydrogenBondsPairsCreator,
                new HbondsPunishHOCAngleCreator(hydrogenBondsCreator),
                new HBondsPunishOHNAngleCreator(hydrogenBondsCreator),
        };


        model.atoms().molecularSystem.createDistanceMatrix();
        System.out.println(" --------------------- Add atoms starts ----------------------------------------") ;
//        Utils.addAtoms(model, false, commands,
//                new PlaneParametersList(EnergyCreator.parametersDirectory(commands) +
//                        "/" + MeshiPotential.PLANE_PARAMETERS), false);
        Utils.addHydrogens(model, commands);
        System.out.println(" --------------------- Add atoms done ----------------------------------------") ;

        DistanceMatrix distanceMatrix = model.atoms().molecularSystem.getDistanceMatrix();
        distanceMatrix.update(1);
        TotalEnergy minimizationEnergy = new TotalEnergy(model,distanceMatrix,energyCreators,commands,model.atoms(),"My energy function.");
        minimizationEnergy.evaluate();

        HydrogenBondsEnergy hydrogenBondsEnergy = (HydrogenBondsEnergy) hydrogenBondsCreator.term();
        System.out.println(hydrogenBondsEnergy);
        System.out.println("\n\n H-Bonds list\n");
        HB_DistanceAttribute attribute;
        for (Distance distance : hydrogenBondsEnergy.hBondList())
                if (accept(distance)) {
                    display(distance);
                    attribute = (HB_DistanceAttribute) distance.getAttribute(HB_DistanceAttribute.key);
                    System.out.println(attribute.isHbond);
                }

        System.out.println("\n\n H-Bond pairs list\n");
        System.out.println("Size of (HydrogenBondsPairsEnergy) hydrogenBondsPairsCreator.term()).getPairsOfHBEElementsList() = "+((HydrogenBondsPairsEnergy) hydrogenBondsPairsCreator.term()).getPairsOfHBEElementsList().size());
        int i = 0;
        for (PairOfHydrogenBondsElements element : ((HydrogenBondsPairsEnergy) hydrogenBondsPairsCreator.term()).getPairsOfHBEElementsList()){
            Distance hb1 = element.HOelement1();
            Distance hb2 = element.HOelement2();
            System.out.println("Distances "+i+" are: \n"+hb1+"\n"+hb2+"\n");
            i++;
            if (accept(hb1) & accept(hb2)) {
                display(hb1);
                display(hb2);
                System.out.println("goodPair: " + element.isGoodPair());
                System.out.println();
            }
        }
        PairsOfHBEElementsList pairsOfHBEElementsList = ((HydrogenBondsPairsEnergy) hydrogenBondsPairsCreator.term()).getPairsOfHBEElementsList();
        String pdbName = model.toString().substring(0, 5);
        patternStringMatrixPairsUpdater(pairsOfHBEElementsList, patternMatrix, stringMatrix, pdbName, pairsLine);
        System.out.println(pairsLine[0]);
        //patterMatrix[0][0]++;
        //exportPairOfHydrogenBondsElements(hydrogenBondsPairsCreator, model);
        //int numberOfHBP = exportGoodPairOfHydrogenBondsElements(hydrogenBondsPairsCreator, model);
        //boolean[][] commonHBMatrix = commonHBMatrixMaker (hydrogenBondsPairsCreator, numberOfHBP, model);

      }


    private static boolean accept(Distance distance) {
         return true;
//         int[] range = {15,28};
//         Atom atom1 = distance.atom1();
//         Atom atom2 = distance.atom2();
//         int n1 = atom1.residueNumber();
//         int n2 = atom2.residueNumber();
//         double bondLength = distance.distance(); /*the bond length*/
//         if ((n1 >= range[0]) && (n2 >= range[0]) && (n1 <= range[1]) && (n2 <= range[1]) /*bond length criteria*/ && (bondLength<=2.5)) return true;
//         return false;
     }
    private static void display(Distance distance) {
        Atom atom1 = distance.atom1();
        Atom atom2 = distance.atom2();
        System.out.println(distance+" "+atom1.residue()+" "+atom1.type()+" ; "+atom2.residue()+" "+atom2.type());
    }

    private static void exportPairOfHydrogenBondsElements (HydrogenBondsPairsCreator hydrogenBondsPairsCreator, Protein model) throws FileNotFoundException, UnsupportedEncodingException {
        String fileName="nativeStructures/txt/PairOfHydrogenBondsElements." +model.toString()+".txt";
        try (PrintWriter writer = new PrintWriter(fileName, "UTF-8")) {
            writer.println(model);
            int i = 1;
            for (PairOfHydrogenBondsElements element : ((HydrogenBondsPairsEnergy) hydrogenBondsPairsCreator.term()).getPairsOfHBEElementsList()) {
                Distance hb1 = element.HOelement1();
                Distance hb2 = element.HOelement2();

                writer.println(i+" h1: "+hb1.atom1().residueNumber() + " o1: " + hb1.atom2().residueNumber() +
                        " h2: " + hb2.atom1().residueNumber() + " o2: " + hb2.atom2().residueNumber() + " " + element.isGoodPair() + " " + element.getPairFactor());
                i++;
            }
            writer.close();
        }
    }

    private static int exportGoodPairOfHydrogenBondsElements (HydrogenBondsPairsCreator hydrogenBondsPairsCreator, Protein model) throws FileNotFoundException, UnsupportedEncodingException {
        String fileName="nativeStructures/txt/GoodPairOfHydrogenBondsElements." +model.toString()+".txt";
        int i = 1;
        try (PrintWriter writer = new PrintWriter(fileName, "UTF-8")) {
            writer.println(model);
            for (PairOfHydrogenBondsElements element : ((HydrogenBondsPairsEnergy) hydrogenBondsPairsCreator.term()).getPairsOfHBEElementsList()) {
                Distance hb1 = element.HOelement1();
                Distance hb2 = element.HOelement2();

                if(element.isGoodPair()){
                    writer.println(i+" h1: "+hb1.atom1().residueNumber() + " o1: " + hb1.atom2().residueNumber() +
                            " h2: " + hb2.atom1().residueNumber() + " o2: " + hb2.atom2().residueNumber() + " " + element.isGoodPair() + " " + element.getPairFactor());
                }
                i++;

            }
            writer.close();
        }
        return i-1;
    }

    private static boolean[][] commonHBMatrixMaker (HydrogenBondsPairsCreator hydrogenBondsPairsCreator, int numberOfHBP, Protein model) throws FileNotFoundException, UnsupportedEncodingException {

        boolean[][] commonHBMatrix= new boolean[numberOfHBP][numberOfHBP];
        int i=0;
        for (PairOfHydrogenBondsElements elementI : ((HydrogenBondsPairsEnergy) hydrogenBondsPairsCreator.term()).getPairsOfHBEElementsList()) {
            Distance hb1i = elementI.HOelement1();
            Distance hb2i = elementI.HOelement2();
            int j = 0;
            for (PairOfHydrogenBondsElements elementJ : ((HydrogenBondsPairsEnergy) hydrogenBondsPairsCreator.term()).getPairsOfHBEElementsList()) {
                Distance hb1j = elementJ.HOelement1();
                Distance hb2j = elementJ.HOelement2();
                if (i == j) {
                    commonHBMatrix[i][j] = false;
                }
                else if (hb1i.equals(hb1j)||hb1i.equals(hb2j)||hb2i.equals(hb1j)||hb2i.equals(hb2j)){
                    commonHBMatrix[i][j] = true;
                }
                j++;
            }
            i++;
        }
        matrixToFile(commonHBMatrix, numberOfHBP, model);
        return commonHBMatrix;
    }

    private static void matrixToFile (boolean[][] commonHBMatrix, int numberOfHBP, Protein model) throws FileNotFoundException, UnsupportedEncodingException {
        String fileName="nativeStructures/txt/commonHBMatrix." +model.toString()+".txt";
        boolean isSymmetric = true;
        try (PrintWriter writer = new PrintWriter(fileName, "UTF-8")) {
            /*writer.print(" ");
            int i=1;
            while(i<=numberOfHBP){
                writer.print(i);
                i++;
            }
            writer.println();*/
            //int j = 1;
            int i=0;
            for (boolean[] elementLine : commonHBMatrix) {
                //writer.print(j);
                int j=0;
                for (boolean element : elementLine){
                    if(element){
                        writer.print("1");
                    }
                    else{
                        writer.print("0");
                    }
                    if (isSymmetric && commonHBMatrix[i][j]!=commonHBMatrix[j][i]){
                        isSymmetric =false;
                    }
                    j++;
                }
                writer.println();
                i++;
                //j++;
            }
            writer.close();
        }
        System.out.println("Is the matrix symmetric : "+ isSymmetric);
    }

    private static void patternStringMatrixPairsUpdater (PairsOfHBEElementsList pairsOfHBEElementsList, int[][] patternMatrix, String[][] stringMatrix, String pdbName, String[] pairsLine)
            throws FileNotFoundException, UnsupportedEncodingException, RuntimeException {

        int i = 0;
        int [][] currentPatternMatrix = new int[750][750];
        for (PairOfHydrogenBondsElements elementI : pairsOfHBEElementsList) {
            Distance hb1i = elementI.HOelement1();
            Distance hb2i = elementI.HOelement2();
            if (elementI.isGoodPair()) {  // added this if statement isGoodPair
                int j = 0;
                for (PairOfHydrogenBondsElements elementJ : pairsOfHBEElementsList) {
                    Distance hb1j = elementJ.HOelement1();
                    Distance hb2j = elementJ.HOelement2();
                    if (i < j && elementJ.isGoodPair() && (hb1i.equals(hb1j) || hb1i.equals(hb2j) || hb2i.equals(hb1j) || hb2i.equals(hb2j))) {  //i<j was i!=j added isGoodPair
                        int matrixI = findPattern(elementI);
                        int matrixJ = findPattern(elementJ);
                        if (matrixI != -1 && matrixJ != -1) {
                            patternMatrix[matrixI][matrixJ]++;
                            currentPatternMatrix[matrixI][matrixJ]++;
                            if (stringMatrix[matrixI][matrixJ] == null) {
                                stringMatrix[matrixI][matrixJ] = "";
                            }
                            //*********************
                            String toWrite=pdbName+".";
                            if (hb1i.atom1.atom.name=="H" && hb1i.atom2.atom.name =="O"){
                                toWrite=toWrite+hb1i.atom1().pdbLine().residueNumber() + "." + hb1i.atom2().pdbLine().residueNumber() + ".";
                            }
                            else if (hb1i.atom2.atom.name=="H" && hb1i.atom1.atom.name =="O"){
                                toWrite=toWrite+hb1i.atom2().pdbLine().residueNumber() + "." + hb1i.atom1().pdbLine().residueNumber() + ".";
                            }
                            else{
                                throw new RuntimeException("Weird atoms in the HBondPairs");
                            }
                            //****
                            if (hb2i.atom1.atom.name=="H" && hb2i.atom2.atom.name =="O"){
                                toWrite=toWrite+hb2i.atom1().pdbLine().residueNumber() + "." + hb2i.atom2().pdbLine().residueNumber() + ".";
                            }
                            else if (hb2i.atom2.atom.name=="H" && hb2i.atom1.atom.name =="O"){
                                toWrite=toWrite+hb2i.atom2().pdbLine().residueNumber() + "." + hb2i.atom1().pdbLine().residueNumber() + ".";
                            }
                            else{
                                throw new RuntimeException("Weird atoms in the HBondPairs");
                            }
                            //****
                            if (hb1j.atom1.atom.name!="H" && hb1j.atom2.atom.name !="O"){
                                toWrite=toWrite+hb1j.atom1().pdbLine().residueNumber() + "." + hb1j.atom2().pdbLine().residueNumber() + ".";
                            }
                            else if (hb1j.atom2.atom.name!="H" && hb1j.atom1.atom.name !="O"){
                                toWrite=toWrite+hb1j.atom2().pdbLine().residueNumber() + "." + hb1j.atom1().pdbLine().residueNumber() + ".";
                            }
                            else {
                                throw new RuntimeException("Weird atoms in the HBondPairs");
                            }
                            //****
                            if (hb2j.atom1.atom.name!="H" && hb2j.atom2.atom.name !="O"){
                                toWrite=toWrite+hb2j.atom1().pdbLine().residueNumber() + "." + hb2j.atom2().pdbLine().residueNumber() + "|";
                            }
                            else if (hb2j.atom2.atom.name!="H" && hb2j.atom1.atom.name !="O"){
                                toWrite=toWrite+hb2j.atom2().pdbLine().residueNumber() + "." + hb2j.atom1().pdbLine().residueNumber() + "|";
                            }
                            else {
                                throw new RuntimeException("Weird atoms in the HBondPairs");
                            }
                            stringMatrix[matrixI][matrixJ] = stringMatrix[matrixI][matrixJ] + toWrite;
                            //***********************
                            /*stringMatrix[matrixI][matrixJ] = stringMatrix[matrixI][matrixJ] + pdbName +
                                    hb1i.atom1().pdbLine().residueNumber() + "." + hb1i.atom2().pdbLine().residueNumber() + "." +
                                    hb2i.atom1().pdbLine().residueNumber() + "." + hb2i.atom2().pdbLine().residueNumber() + "." +
                                    hb1j.atom1().pdbLine().residueNumber() + "." + hb1j.atom2().pdbLine().residueNumber() + "." +
                                    hb2j.atom1().pdbLine().residueNumber() + "." + hb2j.atom2().pdbLine().residueNumber() + "|";*/
                        }
                    }
                    j++;
                }
            }
            i++;
        }
        updateEdgeLine(pairsLine, currentPatternMatrix, pdbName);
    }

    private static void updateEdgeLine (String[] pairsLine, int[][] currentPatternMatrix, String pdbName){
        pairsLine[0] = pdbName+" ";
        for (int i=0; i<currentPatternMatrix.length; i++){
            for (int j=i+1; j<currentPatternMatrix[i].length; j++){
                if (currentPatternMatrix[i][j]>0){
                    // when changing edgeNumber, also change it in graphMaker, graphMakerNoLocations, graphMakerNoLocations
                    int edgeNumber=(i*currentPatternMatrix[i].length)+j;
                    pairsLine[0]=pairsLine[0]+edgeNumber+",";
                }
            }
        }
    }

    private static int findPattern (PairOfHydrogenBondsElements elementI){

        boolean foundI = false;
        int matrixI=-1;

        try {
            BufferedReader br = new BufferedReader(new FileReader("HBBetaHelixParametersSmall.txt"));
            String line = br.readLine();
            int lineNum=0;

            while (!foundI && line != null) {
                if (line.contains("firstBond : "+elementI.seqDistance1()) && line.contains("secondBond : "+elementI.seqDistance2()) &&
                        line.contains("h_dis : "+elementI.hhSeqDistance()) && line.contains("o_dis : "+elementI.ooSeqDistance()) &&
                        line.contains("ho-dis : "+elementI.hoSeqDistance()) && line.contains("oh_dis : "+elementI.ohSeqDistance())){
                    matrixI=lineNum;
                    foundI=true;
                }
                line = br.readLine();
                lineNum++;
            }
            br.close();
        }
        catch (IOException e) {
            e.printStackTrace();
        }

        return matrixI;
    }

    private static void patternMatrixToFile (int[][] patternMatrix) throws FileNotFoundException, UnsupportedEncodingException {
        String fileName="nativeStructures/patternMatrix.txt";

        try (PrintWriter writer = new PrintWriter(fileName, "UTF-8")) {

            int i=0;
            for (int[] elementRow : patternMatrix) {
                int j=0;
                for (int element : elementRow){
                    writer.print(" "+element);
                    j++;
                }
                writer.println();
                i++;
            }
            writer.close();
        }
    }

    private static String[] patternNamesArrayMaker () {

        String[] patternNames = new String[750];
        try {
            BufferedReader br = new BufferedReader(new FileReader("HBBetaHelixParametersSmall.txt"));
            String line = br.readLine();
            int lineNum=0;

            while (line != null) {
                patternNames[lineNum] = patternNameMaker(line, lineNum);
                line = br.readLine();

                lineNum++;
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
        return patternNames;
    }

    private static String patternNameMaker (String line, int lineNum) {
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

    public static String removeSpace (String line) {
        while (line.charAt(0) == ' ') {
            line = line.substring(1);
        }
        return line;
    }

    public static void graphMaker (String[] patternNames, int[][] patternMatrix, String[][] stringMatrix) throws FileNotFoundException, UnsupportedEncodingException {
        String fileName="nativeStructures/graph.txt";

        try (PrintWriter writer = new PrintWriter(fileName, "UTF-8")) {
            writer.println("key HBPattern1 HBPattern2 Weight Locations");
            for (int i=0; i<patternMatrix.length; i++) {
                for (int j=i+1 ; j<patternMatrix[i].length; j++){
                    if (patternMatrix[i][j]>0 && stringMatrix[i][j]!=null){
                        int edgeNumber=(i*patternMatrix[i].length)+j;
                        writer.println(edgeNumber+" "+patternNames[i]+" "+patternNames[j]+" "+patternMatrix[i][j]+" "+stringMatrix[i][j]);
                    }
                }
            }
            writer.close();
        }
    }

    public static void graphMakerNoLocations (String[] patternNames, int[][] patternMatrix, String[][] stringMatrix) throws FileNotFoundException, UnsupportedEncodingException {
        String fileName="nativeStructures/graphNoLocations.txt";

        try (PrintWriter writer = new PrintWriter(fileName, "UTF-8")) {
            writer.println("key HBPattern1 HBPattern2 Weight");
            int i=0;
            for (int[] elementRow : patternMatrix) {
                for (int j=i+1 ; j<elementRow.length; j++){
                    if (elementRow[j]>0 && stringMatrix[i][j]!=null){
                        int edgeNumber=(i*patternMatrix[i].length)+j;
                        writer.println(edgeNumber+" "+patternNames[i]+" "+patternNames[j]+" "+elementRow[j]);
                    }
                }
                i++;
            }
            writer.close();
        }
    }

    public static void locationGraphMaker (int[][] patternMatrix, String[][] stringMatrix) throws FileNotFoundException, UnsupportedEncodingException {
        String fileName="nativeStructures/LocationGraph.txt";

        try (PrintWriter writer = new PrintWriter(fileName, "UTF-8")) {
            writer.println("key Locations");
            int i=0;
            for (int[] elementRow : patternMatrix) {
                for (int j=i+1 ; j<elementRow.length; j++){
                    if (elementRow[j]>0 && stringMatrix[i][j]!=null){
                        int edgeNumber=(i*patternMatrix[i].length)+j;
                        writer.println(edgeNumber+" "+stringMatrix[i][j]);
                    }
                }
                i++;
            }
            writer.close();
        }
    }

    public static void pairsInProteinsMaker (String[] pairsInProteins) throws FileNotFoundException, UnsupportedEncodingException {
        String fileName="nativeStructures/pairsInProteins.txt";

        try (PrintWriter writer = new PrintWriter(fileName, "UTF-8")) {
            writer.println("pdbName pairNumbers");
            int lineNumber=1;
            for (String pairsLine : pairsInProteins) {
                if (pairsLine!=null) {
                    writer.println(pairsLine);
                }
            }
            writer.close();
        }
    }


}
