package programs;

import meshi.energy.EnergyCreator;
import meshi.energy.EvaluationException;
import meshi.energy.TotalEnergy;
import meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma.AtomicPairwisePMFSummaCreator;
import meshi.energy.rg.RgCreator;
import meshi.energy.simpleEnergyTerms.angle.AngleCreator;
import meshi.energy.simpleEnergyTerms.bond.BondCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity.CompositePropensityCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergyCreator;
import meshi.energy.solvate.SolvateCreatorHBforMinimization;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.Atom;
import meshi.sequences.AlignmentException;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.util.*;
import meshi.util.file.MeshiLineReader;
import meshi.util.filters.Filter;
import meshi.util.info.*;

import java.io.File;
import java.io.FileFilter;
import java.io.IOException;
import java.util.ArrayList;

/**

 */
public class SelectModels {
    private static AtomicPairwisePMFSummaCreator summaCreator = new AtomicPairwisePMFSummaCreator();
    private static RamachandranSidechainEnergyCreator ramachCreator = new RamachandranSidechainEnergyCreator();
    private static CompositePropensityCreator propensityCreator = new CompositePropensityCreator();
    private static SolvateCreatorHBforMinimization solvateCreator1 = new SolvateCreatorHBforMinimization(1.0, 0.0, 0.0, 0.0);
    private static SolvateCreatorHBforMinimization solvateCreator2 = new SolvateCreatorHBforMinimization(0.0, 1.0, 0.0, 0.0);
    private static SolvateCreatorHBforMinimization solvateCreator3 = new SolvateCreatorHBforMinimization(0.0, 0.0, 1.0, 0.0);
    private static SolvateCreatorHBforMinimization solvateCreator4 = new SolvateCreatorHBforMinimization(0.01, 0.01, 0.011, 1.0);
    private static RgCreator rgCreator = new RgCreator();
    private static BondCreator bondCreator = new BondCreator();
    private static AngleCreator angleCreator = new AngleCreator();
    private static CompletionFilter completionFilter;

    private static EnergyCreator[] energyCreators = {summaCreator, ramachCreator, propensityCreator,
            solvateCreator1, solvateCreator2, solvateCreator3, solvateCreator4, rgCreator, bondCreator, angleCreator};
    private static QmodElementList qModElements;
    private static File directory;
    private static String key;

    public static final double MIN_RMS = 0.5;

    private static double[] coefficients = null;

    public static void main(String[] argv) throws IOException, UpdateableException,EvaluationException,AlignmentException {
        PdbFilesScanner pdbFilesScanner;
        ProteinInfoList modelInfo;


        CommandList commands = init(argv);

        System.out.println("Analyzing " + directory.getName());
        pdbFilesScanner = new PdbFilesScanner(directory, new OptimizedModelsFilter(), commands);
        modelInfo = pdbFilesScanner.analyze(new AnalyzeDecoys(directory, commands), completionFilter);
        if (coefficients != null) select(coefficients, modelInfo);
        if (modelInfo != null) {
            modelInfo.print();
            //                System.out.println("zScore");
            //                 zScore = modelInfo.toZscore();
            //                 zScore.print();
        }
    }

    private static void select(double[] coefficients, ProteinInfoList modelInfoList) throws AlignmentException{
        double[] gdtij, gdtji;
        double rms;
        DoubleInfoElement iElement;
        int size;
        Protein proteinI, proteinJ;
        for (ProteinInfo info : modelInfoList) {
            if (coefficients.length != info.size() - 1)
                throw new RuntimeException("coefficients.length != info.size()-1");
            double sum = coefficients[0] + coefficients[1] * ((DoubleInfoElement) info.get(2)).value(); //constant + the size-independent qMod
            size = ((IntInfoElement) info.get(0)).value();
            for (int i = 2; i < coefficients.length; i++) {
                iElement = (DoubleInfoElement) info.get(i + 1);
                sum += coefficients[i] * iElement.value() / size;
                info.setSum(sum);
            }
            System.out.println(info.name() + "   " + sum + "    " + name(info));
        }

        System.out.println("Starting sort");
        modelInfoList.sort();
        System.out.println("Sort done");
        for (ProteinInfo info : modelInfoList)
            System.out.println(info.name() + "   " + info.sum() + "    " + name(info));

        modelInfoList.filterByOriginalModel();
        int onCounter = 0;
        for (ProteinInfo proteinInfo : modelInfoList) {
            if (onCounter > 50) proteinInfo.turnOff();
            else if (proteinInfo.on()) onCounter++;
        }
        System.out.println("Selected");
        for (int i = 0; i < modelInfoList.size(); i++) {
            ProteinInfo infoI = modelInfoList.get(i);
            if (infoI.on()) {
                System.out.println(infoI.name() + "   " + infoI.sum() + "    " + name(infoI));
                proteinI = Protein.getCAproteinFromApdbFile(new File(infoI.fileName()));
                for (int j = i + 1; j < modelInfoList.size(); j++) {
                    ProteinInfo infoJ = modelInfoList.get(j);
                    if (infoJ.on()) {
                        proteinJ = Protein.getCAproteinFromApdbFile(new File(infoJ.fileName()));
                        gdtij = Rms.gdt(proteinI, proteinJ);
                        gdtji = Rms.gdt(proteinJ, proteinI);
                        try {
                            rms = Rms.rms(proteinI, proteinJ, ResidueAlignmentMethod.BY_RESIDUE_NUMBER);
                        } catch (Exception ex){throw new RuntimeException(ex.getMessage());}
                            if ((gdtij[0] == 1) & (gdtji[0] == 1)) {
                            if (rms < MIN_RMS) {
                                modelInfoList.get(j).turnOff();
                            }
                        }
                    }
                }
            }
        }
    }

    private static double[] getCoefficients(String fileName) throws IOException {
        MeshiLineReader reader = new MeshiLineReader(fileName);
        String line = reader.readLine();
        String[] words = line.split(",");
        double[] coefficients = new double[words.length];
        int i = 0;
        for (String word : words) {
            coefficients[i] = Double.valueOf(word);
            i++;
        }
        return coefficients;
    }

    private static class AnalyzeDecoys implements ProteinAnalyzer {
        Protein nativeStructure;
        CommandList commands;
        File directory;

        public AnalyzeDecoys(File directory, CommandList commands) {
            File[] temp = directory.listFiles(new NativeFilter());
            if (temp.length == 0) nativeStructure = null;
            else {
                nativeStructure = Protein.getCAproteinFromApdbFile(temp[0]);
            }
            this.commands = commands;
            this.directory = directory;
        }

        public ProteinInfo analyze(Protein model) throws UpdateableException, EvaluationException,AlignmentException {
            File dsspFile = getDsspFile(model.sourceFile());
            Utils.AssignDSSP(model, dsspFile.getAbsolutePath());
            ProteinInfo out = new ProteinInfo(model.name(), model.sourceFile().getAbsolutePath(), "ModelInfoFor_" + directory.getName(), model);
            out.add(new IntInfoElement(InfoType.SIZE, "Number of model atoms", model.atoms().size()));
            double[] gdt = Rms.gdt(nativeStructure, model);
            out.add(new DoubleInfoElement(InfoType.GDT_TS, "GDT_TS from native", gdt[0]));
            out.add(qModElements.getElement(model.name().substring(0, model.name().indexOf('.'))));
            TotalEnergy energy = new TotalEnergy(model, energyCreators, commands,"generic energy (analyze)");
            try {
                energy.evaluate();
            } catch (Exception ex) {
                System.out.println(" Failed to evaluate energy  for " + model + "due to" + ex);
                return null;
            }
            out.add(energy.energyInfo());
            DistanceMatrix distanceMatrix = model.atoms().molecularSystem().getDistanceMatrix();
            out.add(new IntInfoElement(InfoType.NBL, "Non-bonded list size", distanceMatrix.nonBondedList().size()));

            return out;
        }

        private void removeBackboneAtoms(Protein model) {
            for (Atom atom : model.atoms()) {
                if (atom.type().backbone()) atom.resetCoordinates();
            }
        }

        public ProteinInfoList getProteinInfoList() {
            return new ProteinInfoList("ModelInfoFor_" + directory.getName());
        }
    }

    //------------------------------------------- filters --------------------------------------------------------------------------------------------------------

    private static class NativeFilter implements FileFilter {
        public boolean accept(File file) {
            if (file.getName().equals("native.pdb")) return true;
            return false;
        }
    }

    private static class OptimizedModelsFilter implements FileFilter {
        public boolean accept(File file) {
            String name = file.getName();
            String[] nameElements = name.replace('.', '%').split("%");
            int nameElementsLength = nameElements.length;
            if (nameElementsLength < 3) return false;
            if (!nameElements[nameElementsLength - 1].equals("pdb")) return false;
            if (nameElements[nameElementsLength - 2].equals("out")) return true;
            if (nameElements[nameElementsLength - 3].equals("out")) return true;
            return false;
        }
    }

    private static class CaspTargetDirectoryFilter implements FileFilter {
        public boolean accept(File file) {
            String name = file.getName();
            if ((name.length() == 5) & (name.charAt(0) == 'T') & (name.charAt(1) == '0')) return true;
            else return false;
        }
    }

    //------------------------------------------------------------ qMod ----------------------------------------------------

    private static class QmodElement extends DoubleInfoElement {
        String modelName;

        public QmodElement(String line) {
            super(InfoType.QMOD, "model quality estimate", Double.valueOf(line.substring(line.indexOf('\t') + 1)));
            modelName = line.substring(0, line.indexOf('\t'));
        }
    }

    private static class QmodElementList extends ArrayList<QmodElement> {
        public QmodElementList(File file) throws IOException {
            if (file.exists()) {
                MeshiLineReader reader = new MeshiLineReader(file.getAbsolutePath());
                for (int i = 0; i < 4; i++) reader.readLine();
                String line;
                while ((line = reader.readLine()) != null)
                    if (!line.equals("END")) add(new QmodElement(line));
            }
        }

        public QmodElement getElement(String name) {
            for (QmodElement element : this) {
                if (element.modelName.equals(name)) return element;
            }
            return new QmodElement(name + "\t" + 0.0);
        }
    }

    private static class CompletionFilter implements Filter {
        String line;
        String key;

        public CompletionFilter(String key) {
            this.key = key;
        }

        public boolean accept(Object obj) {
            File file = (File) obj;
            try {
                MeshiLineReader reader = new MeshiLineReader(file);
                while ((line = reader.readLine()) != null) {
                    if (line.indexOf(key) != -1) return true;
                }
            } catch (Exception ex) {
                System.out.println("Failed to filter " + file);
                return false;
            }
            return false;
        }
    }

    private static CommandList init(String[] argv) throws IOException {
        if ((argv.length != 4) && (argv.length != 5)) {
            throw new RuntimeException(usageString(" Wrong number of arguments."));
        }
        if (!(new File(argv[0])).exists()) throw new RuntimeException(usageString("Cannot find commands file."));
        if (!(new File(argv[1])).exists()) throw new RuntimeException(usageString("Cannot find models directory"));
        if (!(new File(argv[1])).isDirectory())
            throw new RuntimeException(usageString("Models directory is not a directory"));
        if (!(new File(argv[2])).exists()) throw new RuntimeException(usageString("Cannot find NBT scores"));
        if (argv.length == 5) {
            if (!(new File(argv[4])).exists())
                throw new RuntimeException(usageString("Cannot find regression coefficients."));
            coefficients = getCoefficients(argv[4]);
        }
        CommandList out = new CommandList(argv, new CommandsException("filed to open commands"));
        directory = new File(argv[1]);
        File qmodFile = new File(argv[2]);
        qModElements = new QmodElementList(qmodFile);
        key = argv[3];
        completionFilter = new CompletionFilter(key);

        return out;
    }

    public static File getDsspFile(File file) {
        String path = file.getAbsolutePath();
        int pathLength = path.length();
        String origFileName = path.substring(0, path.indexOf("TS") + 3) + ".pdb";
        //String origFileName1 = path.substring(0,pathLength-8)+".pdb";
        //String origFileName2 = path.substring(0,pathLength-10)+".pdb";
        //String origFileName3 = path.substring(0,pathLength-11)+".pdb";

        File orig = new File(origFileName);
        //if (!orig.exists()) orig = new File(origFileName2);
        //if (!orig.exists()) orig = new File(origFileName3);
        //if (!orig.exists()) throw new RuntimeException("failed to find original file for "+path+" with either "+origFileName1+" , "+origFileName2+" or "+origFileName3);
        if (!orig.exists()) throw new RuntimeException("failed to find original file " + origFileName + " for " + path);
        File dsspFile = new File(orig.getAbsolutePath() + ".dssp");
        if (!dsspFile.exists()) throw new RuntimeException(" failed to find dssp file " + dsspFile.getAbsolutePath());
        return dsspFile;
    }

    private static String usageString(String message) {
        return "Errorr\n" + message + "\n" + "Usage: java -Xmx1000m programs/SelectModels <commandsFileName>  <modelsDirectory> <NBTscores> <complitionKey> [<regressionCoefficients>]";
    }

    private static String name(ProteinInfo info) {
        return ((StringInfoElement) info.get(1)).value();
    }
}
