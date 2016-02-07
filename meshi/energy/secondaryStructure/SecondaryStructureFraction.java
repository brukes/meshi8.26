package meshi.energy.secondaryStructure;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyInfoElement;
import meshi.energy.EnergyType;
import meshi.energy.TotalEnergy;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.parameters.SecondaryStructure;

/**
 * Created by IntelliJ IDEA.
 * User: chen
 * Date: 20/05/12
 * Time: 11:26
 * To change this template use File | Settings | File Templates.
 */
public class SecondaryStructureFraction extends AbstractEnergy {
    Protein protein;
    public SecondaryStructureFraction(EnergyInfoElement info, Protein protein){
        super(toArray(),info, EnergyType.NON_DIFFERENTIAL);
        this.protein = protein;
    }

    public EnergyInfoElement evaluate() {
        int sum = 0;
        double sumSS = 0;
        for (Residue residue : protein.residues()) {
            if (!residue.dummy()) {
                sum++;
                if (residue.secondaryStructure() != SecondaryStructure.COIL)
                    sumSS += 1;
            }
        }
        info().setValue(sumSS/sum);
        return info();
    }

    public void evaluateAtoms() {
    }
    public void test(TotalEnergy te, Atom a) {
             System.out.println("SecondaryStructureFraction tested.");
     }

}

