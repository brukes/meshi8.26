package programs;

import meshi.PDB.PdbReader;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.ResidueList;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.ca.CaResidue;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtoms;
import meshi.sequences.AlignmentException;
import meshi.sequences.ResidueAlignment;
import meshi.sequences.ResidueAlignmentColumn;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.util.GDTcalculator;
import meshi.util.MeshiException;
import meshi.util.Rms;
import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;
import meshi.util.string.StringList;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;

/**

 */
public class PostMortem10 {

    /*
    This is not a real program but rather a set of procedures the needs to be executed
    in a certain order. I simply commented-out the lines that I did not need at a specific stage.
     */
     public static void main(String[] args) throws IOException, AlignmentException {
         /*
            Assumes the working directory is where you want the native structures to be.
          */
         //TargetsPage tp10 = new TargetsPage("casp10");
         //tp10.getNativeStructures();


         /*
         Assumes the working directory includes a sub-directory predictions
         and a sub-directory for each target
          */
         //QA("stage1");
         //QA("stage2");

         //TargetsPage tp10 = new TargetsPage("casp10");
         //tp10.getRefinmentInitialModels();

         //refinementAssessment();

         refinement();

     }

    public static void refinement() throws IOException,AlignmentException{
        MeshiLineReader reader = new MeshiLineReader("servers.txt");
        String line = reader.readLine();
        String[] serverNames = line.split(" ");
        MeshiWriter output = new MeshiWriter("refinement.dat");
        MeshiWriter logger = new MeshiWriter("refinement.log");
        File nativeDir = new File("nativeDomains");
        File[] nativeDomainFiles = nativeDir.listFiles();
        File modelsDir = new File("modelFiles");
        File[] domainDirs = modelsDir.listFiles();

        for (File nativeDomainFile : nativeDomainFiles) {
            AtomList atomList = new AtomList(nativeDomainFile.getAbsolutePath());
            Protein nativeStructure = new Protein(atomList, ResidueExtendedAtoms.creator);
            String domain = nativeDomainFile.getName().substring(0,8);
            File modelDir = null;
            for (File file : domainDirs) {
                if (file.getName().startsWith(domain)) {
                    modelDir = file;
                    break;
                }
            }
            if (modelDir == null) {
                System.out.println("Could not find domain "+domain);
                continue;
            }
            File[] temp = modelDir.listFiles();
            modelDir = temp[0];
            refinementDomain(modelDir,nativeStructure,serverNames,output,logger);
            output.flush();
            logger.flush();
        }
        output.close();
        logger.close();
    }

    public static void refinementDomain(File modelDir,
                                        Protein nativeStructure,String[] serverNames,
                                        MeshiWriter output, MeshiWriter logger)throws AlignmentException{
        File[] files = modelDir.listFiles();
        File[] myModels = new File[5];
        ArrayList<File>serverModels = new ArrayList<File>();
        int modelNumber;

        for (File file : files) {
            if (file.getName().indexOf("315")> 0) {
                modelNumber = getModelNumber(file);
                myModels[modelNumber-1] = file;
            }
            else {
                for (String server : serverNames) {
                    if (file.getName().indexOf(server)>0){
                        serverModels.add(file);
                    }
                }
            }
        }
        refinmentDomainModels(nativeStructure,myModels,serverModels,output,logger);
    }

    public static void refinmentDomainModels(Protein nativeStructure,File[] myModeFiles,
                                             ArrayList<File>serverModelFiles,
                                             MeshiWriter output, MeshiWriter logger) throws AlignmentException {
        double[] origRms = {100000,100000, 10000,10000,10000};
        double rms;
        Protein[] orig = new Protein[5];
        Protein[] myModels = new Protein[5];
        Protein serverModel ;
        AtomList atomList;
        double[] myGDT = new double[5];
        double[] serverGDT = new double[5];
        for (int i = 0; i < 5; i++)  {
            if (myModeFiles[i] != null){
                atomList = new AtomList(myModeFiles[i].getAbsolutePath());
                myModels[i] = new Protein(atomList, ResidueExtendedAtoms.creator);
                if (myModels[i].chain().numberOfNonDummyResidues() < 10) myModels[i] = null;
                    myGDT[i] = Rms.gdt(nativeStructure,myModels[i],
                                GDTcalculator.Type.TS,ResidueAlignmentMethod.BY_RESIDUE_NUMBER)[0];
            }
            else {
                myModels[i] = null;
                myGDT[i] = -1;
            }
        }
        int k = 0;
        for (File serverModelFile : serverModelFiles){
            System.out.print('.');
            if (k%100==0)System.out.println();
            k++;
            atomList = new AtomList(serverModelFile.getAbsolutePath());
            serverModel = new Protein(atomList, ResidueExtendedAtoms.creator);
            int nonDummy = serverModel.chain().numberOfNonDummyResidues();
            if (nonDummy > 10) {
                for (int i = 0; i < 5; i++)  {
                    if (myModels[i] != null) {
                        try {
                            rms = Rms.rms(myModels[i],serverModel,
                                      ResidueAlignmentMethod.BY_RESIDUE_NUMBER);
                            if (rms < origRms[i]) {
                                origRms[i] = rms;
                                orig[i] = serverModel;
                            }
                        } catch (Exception ex) {continue;};

                    }
                }
            }
        }
        for (int i = 0; i < 5; i++)  {
            if (orig[i] != null){
                serverGDT[i] = Rms.gdt(nativeStructure,orig[i],
                                    GDTcalculator.Type.TS,ResidueAlignmentMethod.BY_RESIDUE_NUMBER)[0];
                System.out.println("\n"+myModeFiles[i].getName()+"   "+myGDT[i]+" "+orig[i].name()+" "+
                        origRms[i]+" "+serverGDT[i]+" "+(myGDT[i]-serverGDT[i]));
                logger.println(myModeFiles[i].getName()+"   "+myGDT[i]+" "+orig[i].name()+" "+origRms[i]+
                        " "+serverGDT[i]+" "+(myGDT[i]-serverGDT[i]));
                output.println(serverGDT[i]+" "+" "+origRms[i]+" "+(myGDT[i]-serverGDT[i]));
            }
        }
    }

    public static int getModelNumber(File file) {
        String s = file.getName().substring(20,21);
        return (new Integer(s).intValue());
    }

    public static void FM_refinment() throws IOException,AlignmentException{
        MeshiLineReader reader = new MeshiLineReader("list");
        String line;
        while ((line = reader.readLine())!= null){
            StringTokenizer tokenizer = new StringTokenizer(line);
            String nativeFileName  = tokenizer.nextToken();
            String serverModelName = tokenizer.nextToken();
            String myModelName         = tokenizer.nextToken();
            AtomList atomList = new AtomList(nativeFileName);
            Protein nativeStructure = new Protein(atomList, CaResidue.creator);
            atomList = new AtomList(serverModelName);
            Protein serverModel = new Protein(atomList, CaResidue.creator);
            atomList = new AtomList(myModelName);
            Protein myModel = new Protein(atomList, CaResidue.creator);
            double gdtOld = Rms.gdt(nativeStructure,serverModel,
                    GDTcalculator.Type.TS,ResidueAlignmentMethod.BY_RESIDUE_NUMBER)[0];
            double gdtNew = Rms.gdt(nativeStructure,myModel,
                    GDTcalculator.Type.TS,ResidueAlignmentMethod.BY_RESIDUE_NUMBER)[0];
            double rmsOld, rmsNew,change;
            try {
                rmsOld = Rms.rms(nativeStructure,serverModel,
                                    ResidueAlignmentMethod.BY_RESIDUE_NUMBER);
                rmsNew = Rms.rms(nativeStructure,myModel,
                                    ResidueAlignmentMethod.BY_RESIDUE_NUMBER);
                change = Rms.rms(serverModel,myModel,ResidueAlignmentMethod.BY_RESIDUE_NUMBER);
            } catch (Exception ex) {continue;}
            System.out.printf("%-10s %-50s %-20s %12.3f %12.3f %12.3f %12.3f\n",nativeFileName,serverModelName,myModelName,
                    change,gdtOld, gdtNew-gdtOld,rmsNew-rmsOld);

        }
    }
    public static void refinementAssessment() throws IOException,AlignmentException{
        File here = new File(".");
        File[] dirList = here.listFiles();
        MeshiWriter writer = new MeshiWriter("refinementResults");
        for (File targetDir : dirList) {
            if (targetDir.isDirectory() && targetDir.getName().startsWith("TR")) {
                String targetName = targetDir.getName();
                System.out.println(targetName);
                String nativeFileName = "..\\nativeStructures\\T0"+targetName.substring(2)+".N.pdb";
                File nativeFile = new File(nativeFileName);
                if (nativeFile.exists()) {
                    AtomList atomList = new AtomList(nativeFileName);
                    Protein nativeStructure =  new Protein(atomList, CaResidue.creator);
                    File[] modelFiles = targetDir.listFiles();
                    Protein initialModel = null;
                    ArrayList<Protein> models = new ArrayList<Protein>();
                    for (File modelFile : modelFiles) {
                        System.out.println(modelFile.getAbsolutePath());
                        atomList = new AtomList(modelFile.getAbsolutePath(),true);
                        if(modelFile.getName().endsWith("init.pdb"))
                            initialModel = new Protein(atomList, CaResidue.creator);
                        else {
                            Protein model = new Protein(atomList, CaResidue.creator);
                            models.add(model);
                        }
                    }
                    ResidueAlignment alignment = new ResidueAlignment(nativeStructure.chain(),
                                                                  "native",
                                                                  initialModel.chain(),
                                                                  "initialModel",
                                                                  ResidueAlignmentMethod.IDENTITY);
                    ResidueList targetResidues = new ResidueList();
                    for (ResidueAlignmentColumn column : alignment){
                        targetResidues.add(column.residue0());
                    }

                    atomList = new AtomList(targetResidues.get(0).atoms().get(0));
                    for (Residue residue : targetResidues)
                        for (Atom atom : residue.atoms())
                            atomList.add(atom);
                    Protein targetStructure = new Protein(atomList, CaResidue.creator);

                    writer.print(targetName+" "+Rms.gdt(targetStructure,initialModel, GDTcalculator.Type.TS,ResidueAlignmentMethod.IDENTITY)[0]);
                    for (Protein model : models)
                        writer.print("  "+ Rms.gdt(targetStructure,model,GDTcalculator.Type.TS,ResidueAlignmentMethod.IDENTITY)[0]);
                    writer.println();
                }
            }
        }
        writer.close();
    }

    public static void QA(String stage) throws IOException,AlignmentException{
        File nativeStructuresDir = new File("..\\..\\nativeStructures");
        File[] nativeStructures = nativeStructuresDir.listFiles();
        File predictionsDir = new File("predictions");
        File[] predictionFiles = predictionsDir.listFiles();
        for (File nativeStructureFile : nativeStructures) {
            Protein nativeStructure = new Protein(new AtomList(new PdbReader(nativeStructureFile.getAbsolutePath())),
                    CaResidue.creator);
            String fullName = nativeStructureFile.getName();
            String name = fullName.substring(0,5);
            System.out.println(name);
            File modelDir = new File(name+"/serverModels");
            if (modelDir.exists()) {
            File[] modelFiles = modelDir.listFiles();
            for (File prediction : predictionFiles) {
                if (prediction.getName().startsWith(name)) {
                    File outputFile   = new File(name+"."+stage+".postMortem");
                    MeshiWriter writer;
                    if ((!outputFile.exists()) & (outputFile.length()==0)) {
                        writer = new MeshiWriter(outputFile);
                        StringList lines = new StringList(prediction);
                        for (String line : lines) {
                            String[] words = line.split("\t");
                            if ((words.length > 1) &&(!words[1].equals("X"))) {
                                for (File modelFile : modelFiles) {
                                    String fileName = modelFile.getName();
                                    if(fileName.indexOf(words[0])>= 0) {
                                        Protein model = new Protein(new AtomList(new PdbReader(modelFile.getAbsolutePath())),
                                                CaResidue.creator);
                                        try {
                                            double[] gdt = Rms.gdt(model, nativeStructure,GDTcalculator.Type.TS, ResidueAlignmentMethod.IDENTITY);
                                            //System.out.println(model.name()+"\t"+gdt[0]+" "+words[1]);
                                            writer.println(model.name()+"\t"+gdt[0]+" "+words[1]);
                                        }
                                        catch (MeshiException me){
                                            System.out.println("Cannot calculate gdt of "+model+" and "+nativeStructure);
                                        }
                                    }
                                }
                            }
                        }
                        writer.close();
                    }
                }
            }
            }
            else
                System.out.println(modelDir.getName()+" does not exist");
        }
    }
}
