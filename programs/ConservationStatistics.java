package programs;

import meshi.energy.EnergyInfoElement;
import meshi.energy.conservationContacts.ConservationContactsCreator;
import meshi.energy.conservationContacts.ConservationContactsRatio;
import meshi.molecularElements.ConSeq;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.sequences.AlignmentException;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.util.CommandList;
import meshi.util.CommandsException;
import meshi.util.Rms;
import meshi.util.info.InfoType;

import java.io.File;
import java.io.IOException;

/**

 */
public class ConservationStatistics {
    private static Protein nativeStructure = null;
    public static void main(String[] argv) throws IOException,AlignmentException {
        File  dir = new File(argv[0]);
        File[] files = dir.listFiles();

        for (File file: files) {
            String fileName = file.getAbsolutePath();
            if (fileName.endsWith("pdb")) {
                System.out.println(fileName);
                Protein protein                     =   getProtein(file);
                getNative(file);
                ConservationContactsRatio score = getConservationContactsRatio(protein,fileName);
                 double rms = getRms(nativeStructure,protein);
                 double[] gdt = getGdt(nativeStructure,protein);
                EnergyInfoElement info  = score.evaluate();
                double contactRation = info.energy();
                double rgRatio = info.infoList().get(0).value();
                System.out.printf("%-30f %-30f %-30f %-30f\n", contactRation,rgRatio, rms,gdt[0]);
            }
        }
    }
    public static Protein getProtein(File file) {
                new MolecularSystem() ;// Resets atomOne numbers
                AtomList atomList = new AtomList(file.getAbsolutePath());
                new MolecularSystem();
                Protein protein = new Protein(atomList, ResidueExtendedAtomsCreator.creator);
                return protein;
    }
    public static ConservationContactsRatio getConservationContactsRatio(Protein protein, String fileName) throws IOException{
         ConservationContactsCreator creator;
        CommandList commands = new CommandList("commands.txt",new CommandsException("Problem with commands"));
        ConSeq.setConSeq(protein.chain(), fileName); //Assign conservation to the residues
        creator = new ConservationContactsCreator(InfoType.CONTACTS11);
        creator.getWeight(commands);
        creator.setFileFound(true);
        ConservationContactsRatio score = ( ConservationContactsRatio)    creator.createEnergyTerm(protein,null,null);
        return score;
    }

    public static void getNative(File file) {
        String nativeFileName;
        String fileName = file.getAbsolutePath();
        if (fileName.endsWith("N.pdb")) nativeFileName = fileName;
        else {
            String[] subs = fileName.split("\\.");
            int index = fileName.indexOf(subs[subs.length-2]);
            nativeFileName = fileName.substring(0,index)+"N.pdb";
            System.out.println(nativeFileName);
        }
        if ((nativeStructure == null)||
             nativeFileName.indexOf(nativeStructure.name())== -1)
        nativeStructure = getProtein(new File(nativeFileName));
    }


    public static double getRms(Protein protein1, Protein protein2) {
        try {
            return Rms.rms(protein1,protein2,ResidueAlignmentMethod.IDENTITY, Rms.RmsType.CA);// protein1,protein2,alignment);
        }catch (Exception ex ) {throw new RuntimeException(ex.getMessage());}
    }
    public static double[] getGdt(Protein protein1, Protein protein2)throws AlignmentException {
        return Rms.gdt(protein1,protein2);
    }

}
