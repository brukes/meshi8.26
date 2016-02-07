package programs;

import meshi.molecularElements.Residue;
import meshi.parameters.SecondaryStructure;
import meshi.sequences.AlignmentException;
import meshi.util.Rms;
import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;
import meshi.util.filters.Filter;
import meshi.util.Utils;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.energy.rg.*;
import meshi.energy.rg.filters.*;

import java.io.IOException;
import java.util.Arrays;

/**
 */
public class RgStataistics {
    public static void main(String[] args) throws IOException, AlignmentException{
        MeshiLineReader reader = new MeshiLineReader(args[0]);
        MeshiWriter writer = new MeshiWriter(args[1]);
        MeshiWriter RgStatistics = new MeshiWriter("RgSatistics");
        MeshiWriter dictionaryWriter = new MeshiWriter("Dictionary.txt");
        String fileName;


        Filter hydrophobicSideChains = new HydrophobicSideChains();
        Filter backboneFilter = new BackboneCAs();
        Filter polarFilter = new PolarSideChains();
        Filter heavyAtomsFilter = new HeavyAtoms();
        Filter secondaryStructureFilter = new SecondaryStructureFilter();
        Filter coilFilter = new CoilFilter();
        Protein protein;
        Runtime runtime = Runtime.getRuntime();
        int i = 0;
        String nativeFileName = null;
        String oldNativeFileName = null;
        Protein nativeStructure = null;
        double[] gdt = {0,0,0,0,0};
        MeshiWriter targetRG = null,targetNames = null;
        while ((fileName = reader.readLine()) != null) {
            if(fileName.length() > 20) {// that is a decoy
                nativeFileName = fileName.substring(0,36)+"N.pdb";
                if (!nativeFileName.equals(oldNativeFileName)) {
                    if (targetNames != null) targetNames.close();
                    if (targetRG != null) targetRG.close();
                    oldNativeFileName = nativeFileName;
                    nativeStructure = new Protein(new AtomList(nativeFileName), new ResidueExtendedAtomsCreator());
                    targetNames = new MeshiWriter(fileName.substring(28,35)+"_names.dat");
                    targetRG = new MeshiWriter(fileName.substring(28,35)+"_RG.dat");
                }
            }


            new MolecularSystem();
            protein = new Protein(new AtomList(fileName), new ResidueExtendedAtomsCreator());
            if (nativeStructure != null)
                gdt = Rms.gdt(nativeStructure,protein);
            removeOutliers(protein);


            String dsspFileName;
            if (fileName.length() < 20)
                dsspFileName = "dssp/" + protein.name() + ".dssp";
            else dsspFileName = "dssp/" + protein.name().substring(0,8)+"/"+protein.name() + ".dssp";
            if (! Utils.AssignDSSP(protein, dsspFileName)) continue;
            for (int iResidue = 0; iResidue< protein.residues().size(); iResidue++) {
                Residue residue = protein.residues().get(iResidue);
                if (residue.secondaryStructure() == SecondaryStructure.UNK) {
                        residue.setSecondaryStructure(SecondaryStructure.COIL);
                            }
            }
            AtomList heavyAtoms = protein.chain().atoms().filter(heavyAtomsFilter);
            System.out.println("analysing " + protein.name());
            AtomList nonPolarSS = protein.chain().atoms().filter(hydrophobicSideChains).filter(secondaryStructureFilter);
            AtomList nonPolarCoil = protein.chain().atoms().filter(hydrophobicSideChains).filter(coilFilter);
            if (nonPolarSS.size() < 50) continue;
            if (nonPolarCoil.size() < 10) continue;
            AtomList backboneSS = protein.chain().atoms().filter(backboneFilter).filter(secondaryStructureFilter);
            AtomList backboneCoil = protein.chain().atoms().filter(backboneFilter).filter(coilFilter);
            if (backboneSS.size() < 10) continue;
            if (backboneCoil.size() == 10) continue;
            AtomList polarSS = protein.chain().atoms().filter(polarFilter).filter(secondaryStructureFilter);
            AtomList polarCoil = protein.chain().atoms().filter(polarFilter).filter(coilFilter);
            if (polarSS.size() < 10) continue;
            if (polarCoil.size() < 10) continue;
            //
            //
            double[] hydrophobicSideChainsSsMahalanobis = getMahalanobis(nonPolarSS);
            double[] hydrophobicSideChainsSsRG = getRG(nonPolarCoil);
            double[] polarsSsMahalanobis = getMahalanobis(polarSS);
            double[] polarsSsRG = getRG(polarCoil);
            double[] hydrophobicSideChainsSsMContacts8 = getContacts(nonPolarSS,8);
            double[] polarSsMContacts8 = getContacts(polarSS,8);
            double[] hydrophobicSideChainsSsMContacts6 = getContacts(nonPolarSS,6);
            double[] polarSsMContacts6 = getContacts(polarSS,6);
            double[] hydrophobicSideChainsSsMContacts10 = getContacts(nonPolarSS,10);
            double[] polarSsMContacts10 = getContacts(polarSS,10);
            double[] hydrophobicSideChainsSsMContacts12 = getContacts(nonPolarSS,12);
            double[] polarSsMContacts12 = getContacts(polarSS,12);
            double[] hydrophobicSideChainsSsMContacts14 = getContacts(nonPolarSS,14);
            double[] polarSsMContacts14 = getContacts(polarSS,14);
            double fragmentedSS = getFragmentedSS(protein.chain().atoms().filter(secondaryStructureFilter));

            dictionaryWriter.println(i + "\t" + fileName);
            if (targetNames != null) targetNames.println(i + "\t" + fileName);
            String out = i + "\t" +   //1
                    gdt[0]+"\t"+ //2
                    Math.log(nonPolarSS.size()) + "\t" + //3
                    Math.log(nonPolarCoil.size()) + "\t" + //4
                    Math.log(polarSS.size()) + "\t" + //5
                    Math.log(polarCoil.size()) + "\t" + //6
                    hydrophobicSideChainsSsMContacts6[0]+"\t"+      //7
                    hydrophobicSideChainsSsMContacts6[1]+"\t"+      //8
                    polarSsMContacts6[0]+"\t"+ //9
                    polarSsMContacts6[1]+"\t"+ //10
                    hydrophobicSideChainsSsMContacts8[0]+"\t"+ //11
                    hydrophobicSideChainsSsMContacts8[1]+"\t"+ //12
                    polarSsMContacts8[0]+"\t"+  //13
                    polarSsMContacts8[1]+"\t"+  //14
                    hydrophobicSideChainsSsMContacts10[0]+"\t"+ //15
                    hydrophobicSideChainsSsMContacts10[1]+"\t"+ //16
                    polarSsMContacts10[0]+"\t"+  //17
                    polarSsMContacts10[1]+"\t"+ // 18
                    hydrophobicSideChainsSsMContacts12[0]+"\t"+ //19
                    hydrophobicSideChainsSsMContacts12[1]+"\t"+ //20
                    polarSsMContacts12[0]+"\t"+  //21
                    polarSsMContacts12[1]+"\t"+ // 22
                    hydrophobicSideChainsSsMContacts14[0]+"\t"+ //23
                    hydrophobicSideChainsSsMContacts14[1]+"\t"+ //24
                    polarSsMContacts14[0]+"\t"+  //25
                    polarSsMContacts14[1]+"\t"+ // 26
                    fragmentedSS; //27


            writer.println(out);
            if (targetRG != null) targetRG.println(out);

            writer.flush();
            dictionaryWriter.flush();
            runtime.gc();
            i++;
        }
        writer.close();
        dictionaryWriter.close();

        if (targetNames != null) targetNames.close();
        if (targetRG != null) targetRG.close();
    }

    private static double[] getRG(AtomList atomList) {
        CenterOfMass cm = new CenterOfMass(atomList);
        return getRG(atomList,cm);
    }

    private static double[] getRG(AtomList atomList, CenterOfMass cm) {
        double[] out = new double[2];
        double[] distances = new double[atomList.size()];
        double d2, dx, dy, dz, rg = 0, minD = 100000;
        int i = 0;
        for (Atom atom : atomList) {
            dx = atom.x() - cm.getX();
            dy = atom.y() - cm.getY();
            dz = atom.z() - cm.getZ();
            d2 = dx * dx + dy * dy + dz * dz;
            distances[i] = d2;
            i++;
        }
        Arrays.sort(distances);
        double sum = 0;
        for (int j = 0; j < atomList.size()/2; j++)
            sum += distances[j];
        out[0] = Math.sqrt(sum / atomList.size()/2);
        for (int j = atomList.size(); j < atomList.size(); j++)
            sum += distances[j];
        out[1] = Math.sqrt(sum / atomList.size());
        return out;
    }

    private static double getAD(AtomList atomList, CenterOfMass cm) {
        double dx, dy, dz, rg = 0;
        for (Atom atom : atomList) {
            dx = atom.x() - cm.getX();
            dy = atom.y() - cm.getY();
            dz = atom.z() - cm.getZ();
            rg += Math.sqrt(dx * dx + dy * dy + dz * dz);
        }
        return rg / atomList.size();
    }

    public static void removeOutliers(Protein protein) {
        for (Atom atom:protein.atoms())  {
            if ((!atom.nowhere()) && (atom.x() < -999))
                atom.setNowhere();
        }
    }

    public static double[] getContacts(AtomList atoms, double threshold) {
        int[] nContacts = new int[atoms.size()];
        for (int i = 0; i < atoms.size(); i++)  {
            Atom atomI = atoms.get(i);
            for (int j = i+1; j < atoms.size();j++) {
                Atom atomJ = atoms.get(j);
                if (atomI.distanceFrom(atomJ) <= threshold) {
                    nContacts[i]--;
                    nContacts[j]--;
                }
            }
        }
        double[] out = new double[2];
        double sum = 0;
        Arrays.sort(nContacts);

        for (int i = 0; i < atoms.size()/2; i++)
            sum += nContacts[i];
        out[0] = sum/(atoms.size()/2);
        for (int i = atoms.size()/2; i < atoms.size()/2; i++)
            sum += nContacts[i];
        out[1] = sum/(atoms.size());

        return out;
    }

    public static double[] getMahalanobis(AtomList nonPolarSS) {
        double[][] coordinateMatrix = new double[nonPolarSS.size()][3];

        int i = 0;
        for (Atom atom : nonPolarSS) {
            coordinateMatrix[i][0] = atom.x();
            coordinateMatrix[i][1] = atom.y();
            coordinateMatrix[i][2] = atom.z();
            i++;
        }
        double[] mu = meshi.util.mathTools.Utils.getMeanVector(coordinateMatrix);
        double[][] covarianceMatrix = meshi.util.mathTools.Utils.getCovarianceMatrix(coordinateMatrix);
        double sum = 0;
        double min = 1000;
        double[] mahal = new double[coordinateMatrix.length];
        double[] out = new double[2];
        for (i = 0; i< coordinateMatrix.length; i++) {
            mahal[i] = Math.sqrt(meshi.util.mathTools.Utils.Mahalanobis(mu, covarianceMatrix, coordinateMatrix[i]));
        }
        Arrays.sort(mahal);
        for (i = 0; i< coordinateMatrix.length/2; i++) {
            sum += mahal[i];
        }
        out[0] = sum/(coordinateMatrix.length/2);
        for (i = coordinateMatrix.length/2; i< coordinateMatrix.length; i++) {
            sum += mahal[i];
        }
        out[1] = sum/(coordinateMatrix.length);
        //if (1 == 1) throw new RuntimeException("???????????????????????") ;
        return out;
    }

    public static double getFragmentedSS(AtomList atoms){
        double out = 0;
        int nResidues = 0;
        SecondaryStructure firstResidueSS = null;
        Residue residue = null;
        int fragmentLength = 0;

        for (Atom atom :atoms){
            if (atom.residue() == residue) continue;
            nResidues++;
            residue = atom.residue();
            if (residue.secondaryStructure() == firstResidueSS) {
                if (firstResidueSS == SecondaryStructure.COIL) continue;
                fragmentLength++;
            }
            else {
                if ((firstResidueSS == SecondaryStructure.HELIX) &&
                        (fragmentLength <= 4)) out += fragmentLength;
                if ((firstResidueSS == SecondaryStructure.SHEET) &&
                        (fragmentLength <= 2)) out += fragmentLength;
                firstResidueSS = residue.secondaryStructure();
                fragmentLength = 1;
            }
        }
        return out/(nResidues+0.00001);
    }
}
