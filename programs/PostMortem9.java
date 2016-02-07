package programs;

import meshi.PDB.PdbLine;
import meshi.energy.EvaluationException;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.sequences.AlignmentException;
import meshi.sequences.ResidueAlignment;
import meshi.sequences.ResidueAlignmentColumn;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.util.*;
import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;
import meshi.util.filters.Filter;
import meshi.util.info.*;
import meshi.util.string.StringList;

import java.io.File;
import java.io.FileFilter;
import java.io.IOException;

/**

 */
public class PostMortem9 {
    private static File rootDirectory;
    private static MeshiWriter serverWriter;
    private static MeshiWriter ourWriter;
    private static String target = "UNKNOWN_TARGET";
    private static ProteinInfoList allOurModels = new ProteinInfoList("allOurModels");
    private static MeshiWriter allModelsWriter;
    private static boolean allModelsWriterFirstTime = true;
    private static double[] zeros = {0.0, 0.0, 0.0, 0.0, 0.0};


    public static void main(String[] argv) throws IOException, UpdateableException, EvaluationException,AlignmentException {

        init(argv);

        System.out.println("Analyzing " + rootDirectory.getName());
        File[] directories = rootDirectory.listFiles();
        for (File directory : directories) {
            if (directory.isDirectory() & (directory.getName().startsWith("T0") || (directory.getName().startsWith("TR")))) {
                System.out.println("Analyzing " + directory);
                analyzeTarget(directory);
            }
        }
    }

    private static void analyzeTarget(File directory) throws IOException, UpdateableException, EvaluationException,AlignmentException {
        File[] files = directory.listFiles();
        if (directoryDone(files)) System.out.println(directory.getName() + " done");
        else {
            serverWriter = new MeshiWriter(directory.getAbsolutePath() + "\\ServerInfo.txt");
            ourWriter = new MeshiWriter(directory.getAbsolutePath() + "\\OurInfo.txt");
            String targetName = directory.getName();
            Protein nativeStructure = getNative(files);
            Protein[] ourModels = getOurModels(files);
            File serverDirectory = getServerDirectory(files);
            PdbFilesScanner serverModels = new PdbFilesScanner(serverDirectory, new ServerModelFilter());
            ProteinInfoList serverModelsInfo = serverModels.analyze(new ServerModelsAnalyzer(nativeStructure, targetName));
            System.out.println(" serverModelsInfo size = " + serverModelsInfo.size());
            Zparameters gdtZparameters = new Zparameters(serverModelsInfo, 1);
            Zparameters rmsZparameters = new Zparameters(serverModelsInfo, 2);
            serverModelsInfo.print(serverWriter);
            analyzeOurModels(ourModels, serverModelsInfo, nativeStructure, gdtZparameters, rmsZparameters, serverDirectory);
            serverWriter.close();
            ourWriter.close();
            (new File(directory.getAbsolutePath() + "\\done")).createNewFile();
        }
    }

    private static boolean directoryDone(File[] files) {
        for (File file : files)
            if (file.getName().equals("done")) return true;
        return false;
    }

    private static class Zparameters {
        private double sum, sum2;
        private double mean, std;
        private int counter;

        public double mean() {
            return sum / counter;
        }

        public double std() {
            return Math.sqrt(sum2 / counter - (sum / counter) * (sum / counter));
        }

        public Zparameters(ProteinInfoList proteinInfoList, int field) {
            sum = 0;
            sum2 = 0;
            counter = 0;
            double value;
            for (ProteinInfo proteinInfo : proteinInfoList) {
                value = proteinInfo.get(field).doubleValue();
                add(value);
            }
            double mean = mean();
            double std = std();
            for (ProteinInfo proteinInfo : proteinInfoList) {
                MeshiInfoElement element = proteinInfo.get(field);
                value = element.doubleValue();
                proteinInfo.add(new DoubleInfoElement("Z_" + element, "Z-score of " + element, (value - mean) / std));
            }


        }

        public void add(double addMe) {
            if (addMe < 9999) {
                sum += addMe;
                sum2 += addMe * addMe;
                counter++;
            }
        }

    }

    private static Protein getNative(File[] files) {
        for (File file : files) {
            if ((file.getName().indexOf("native") != -1) && (file.getName().endsWith(".pdb"))) {
                new MolecularSystem();
                return Protein.getCAproteinFromApdbFile(file);
            }
        }
        throw new RuntimeException(" No native structure has been found in " + files[0].getParent());
    }

    private static Protein[] getOurModels(File[] files) throws IOException {
        Protein[] out = new Protein[5];
        StringList[] stringLists;
        for (File file : files) {
            if (file.getName().endsWith(".casp9")) {
                System.out.println("Extracting models from " + file.getName());
                stringLists = breakToStringLists(file);
                for (int i = 0; i < 5; i++) {
                    //System.out.println("Extracting model "+(i+1)+" using stringList  of size "+stringLists[i].size());
                    out[i] = getProtein(stringLists[i]);
                }
                return out;
            }
        }
        throw new RuntimeException(" Our models are not found in " + files[0].getParent());
    }

    private static Protein getProtein(StringList stringList) {
        MolecularSystem molecularSystem = new MolecularSystem();
        AtomList atoms = new AtomList(molecularSystem);
        String originalModel = "XXXXX1", modelNumber = "XXXXX3";
        for (String string : stringList) {
            if (string.indexOf("attempt") != -1) originalModel = (string.split(" "))[13];
            if (string.indexOf("TARGET") != -1) target = (string.split(" "))[1];
            if (string.indexOf("MODEL") != -1) modelNumber = (string.split(" "))[1];
            if (string.indexOf("ATOM") != -1) {
                atoms.add(new Atom(new PdbLine(string),molecularSystem));
            }
        }
        new MolecularSystem();
        Protein protein = new Protein(atoms, ResidueExtendedAtomsCreator.creator);
        protein.setName(target + "_keasar_" + modelNumber + "_basedOn_" + originalModel);
        return protein;
    }

    private static StringList[] breakToStringLists(File file) throws IOException {
        StringList[] out = new StringList[5];
        MeshiLineReader reader = new MeshiLineReader(file);
        System.out.println(reader.name() + "\t" + reader.readLine());
        String line;
        for (int i = 0; i < 5; i++) {
            out[i] = new StringList();
            boolean done = false;
            while (!done) {
                line = reader.readLine();
                out[i].add(line);
                //System.out.println(line);
                if (line.trim().equals("TER")) done = true;
            }
            System.out.println("out  " + i + " = " + out[i].size());
        }
        reader.close();
        return out;
    }

    private static File getServerDirectory(File[] files) {
        for (File file : files) {
            if (file.isDirectory() && (file.getName().length() == 5) && (file.getName().startsWith("T0") || file.getName().startsWith("TR"))) {
                return file;
            }
        }
        throw new RuntimeException("There is no server models directory in " + files[0].getParent());
    }

    private static class ServerModelFilter implements FileFilter {
        public boolean accept(File file) {
            if (file.getName().endsWith("_TS1")) return true;
            if (file.getName().endsWith("_TS2")) return true;
            if (file.getName().endsWith("_TS3")) return true;
            if (file.getName().endsWith("_TS4")) return true;
            if (file.getName().endsWith("_TS5")) return true;
            if (file.getName().endsWith("_TSX")) return true;
            return false;
        }
    }

    private static class ServerModelsAnalyzer implements ProteinAnalyzer {
        Protein nativeStructure;
        String targetName;

        public ServerModelsAnalyzer(Protein nativeStructure, String targetName) {
            this.nativeStructure = nativeStructure;
            this.targetName = targetName;
        }

        public ProteinInfo analyze(Protein model) throws AlignmentException{
            ProteinInfo out = new ProteinInfo(model.name(), model.sourceFile().getAbsolutePath(), "ModelInfoFor_" + targetName, model);
            out.add(new IntInfoElement(InfoType.SIZE, "Number of model atoms", model.atoms().size()));
            ResidueAlignment residueAignment = new ResidueAlignment(nativeStructure.chain(),
                    nativeStructure.name(),
                    model.chain(), model.name(),
                    ResidueAlignmentMethod.BY_RESIDUE_NUMBER);
            double[] gdt;
            double rms;
            if (residueAignment.size() < 5) {// probably a domain that was not modeled)
                return null;
            } else {
                try  {
                    gdt = Rms.gdt(nativeStructure, model);
                    rms = Rms.rms(nativeStructure, model, ResidueAlignmentMethod.BY_RESIDUE_NUMBER);
                    out.add(new DoubleInfoElement(InfoType.GDT_TS, "GDT_TS from native", gdt[0]));
                    out.add(new DoubleInfoElement(InfoType.RMS, "RMS from native", rms));
                    return out;
                } catch (Exception ex ) {throw new RuntimeException(ex.getMessage());}

            }
        }


        public ProteinInfoList getProteinInfoList() {
            return new ProteinInfoList("ModelInfoFor_" + targetName);
        }
    }

    private static void analyzeOurModels(Protein[] ourModels, ProteinInfoList serverModelsInfo,
                                         Protein nativeStructure, Zparameters gdtZparameters,
                                         Zparameters rmsZparameters,
                                         File serverDirectory) throws IOException,AlignmentException {

        double gdt[], gdtZscore, originalGdt;
        double gdtMean = gdtZparameters.mean();
        double gdtStd = gdtZparameters.std();
        double rms, rmsZscore, originalRms;
        double rmsMean = rmsZparameters.mean();
        double rmsStd = rmsZparameters.std();
        Protein model;
        ProteinInfo serverModelInfo;
        ProteinInfoList modelsInfoList = new ProteinInfoList("OurModels");
        String outputString;
        for (int i = 0; i < 5; i++) {
            model = ourModels[i];
            gdt = Rms.gdt(nativeStructure, model);
            gdtZscore = (gdt[0] - gdtMean) / gdtStd;
            try {
                rms = Rms.rms(nativeStructure, model, ResidueAlignmentMethod.BY_RESIDUE_NUMBER);
            } catch (Exception ex ) {throw new RuntimeException(ex.getMessage());}
            rmsZscore = (rms - rmsMean) / rmsStd;
            serverModelInfo = findOriginalModel(model, serverModelsInfo);
            originalGdt = serverModelInfo.get(1).doubleValue();
            originalRms = serverModelInfo.get(2).doubleValue();
            ProteinInfo modelInfo = new ProteinInfo(model.name(), "target.casp9", model.name(), model);
            modelInfo.add(new DoubleInfoElement(InfoType.GDT_TS, "GDT_TS from native", gdt[0]));
            modelInfo.add(new DoubleInfoElement(InfoType.Z_GDT, "Z-score of GDT_TS", gdtZscore));
            modelInfo.add(new DoubleInfoElement(InfoType.DELTA_GDT_TS, "<GDT_TS of our model> - <GDT_TS of original model> the larger the better", gdt[0] - originalGdt));
            modelInfo.add(new DoubleInfoElement(InfoType.RMS, "RMS from native", rms));
            modelInfo.add(new DoubleInfoElement(InfoType.Z_RMS, "Z-score of RMS", rmsZscore));
            modelInfo.add(new DoubleInfoElement(InfoType.DELTA_RMS, "<RMS of our model> - <RMS of original model> the lower the better", rms - originalRms));
            Protein serverModel = Protein.getExtendedAtomsProteinFromPdbFile(new File(serverDirectory.getAbsolutePath() + "\\" + serverModelInfo.name()));
            modelInfo.add(new DoubleInfoElement("CHANGE",
                    "Structural change during optimization",
                    Rms.rms(model,
                            serverModel,
                            ResidueAlignmentMethod.BY_RESIDUE_NUMBER,
                            new ResidueNumberFilter(model, nativeStructure))));
            modelsInfoList.add(modelInfo);

            if (allModelsWriterFirstTime) {
                allModelsWriter.println(modelInfo.header());
                allModelsWriterFirstTime = false;
            }
            allModelsWriter.println(modelInfo.values());
            allModelsWriter.flush();
            String fileName = serverDirectory.getParent() + "\\" + serverModelInfo.name() + "_keasar_model_" + (i + 1) + ".pdb";
            try {
                model.printAtomsToFile(fileName);
            } catch (Exception ex) {
                throw new RuntimeException("Failed to save file: " + fileName + "due to\n" + ex);
            }
        }
        modelsInfoList.print();
        modelsInfoList.print(ourWriter);
    }

    private static class ResidueNumberFilter implements Filter {
        private int[] residueNumbers;

        public ResidueNumberFilter(Protein protein1, Protein protein2)throws AlignmentException {
            ResidueAlignment residueAlignment = new ResidueAlignment(protein1.chain(), protein1.chain().name(),
                    protein2.chain(), protein2.chain().name(),
                    ResidueAlignmentMethod.BY_RESIDUE_NUMBER);
            residueNumbers = new int[residueAlignment.size()];
            int i = 0;
            for (ResidueAlignmentColumn column : residueAlignment) {
                residueNumbers[i] = column.residue0().number();
                i++;
            }
        }

        public boolean accept(Object obj) {
            ResidueAlignmentColumn column = (ResidueAlignmentColumn) obj;
            int number = column.residue0().number();
            for (int n : residueNumbers) {
                if (n == number) return true;
            }
            return false;
        }
    }


    private static ProteinInfo findOriginalModel(Protein model, ProteinInfoList serverModelsInfo) {
        String name = getOriginalModelName(model);
        for (ProteinInfo serverModelInfo : serverModelsInfo) {
            if (name.equals(serverModelInfo.name())) return serverModelInfo;
        }
        throw new RuntimeException("Failed to find the original model for " + model.name() + " using " + name);
    }

    private static String getOriginalModelName(Protein model) {
        return (model.name().split("_basedOn_"))[1];
    }


    private static void init(String[] argv) throws IOException {
        rootDirectory = new File(".");
        if ((argv.length == 1) && (argv[0].equals("restart"))) {
            resetDirectory();
        }
        allModelsWriter = new MeshiWriter("allOurModels");
    }

    private static void resetDirectory() {
        File[] directories = rootDirectory.listFiles();
        File[] doneFiles = new File[directories.length];
        int i = 0;
        for (File directory : directories) {
            doneFiles[i] = new File(directory.getAbsolutePath() + "\\done");
            i++;
        }
        for (File file : doneFiles) {
            if (file.exists()) file.delete();
        }
    }

    private static String usageString(String message) {
        return "Errorr\n" + message + "\n" + "Usage: java -Xmx1000m programs/SelectModels <commandsFileName>  <modelsDirectory> <NBTscores> <complitionKey> [<regressionCoefficients>]";
    }


}


