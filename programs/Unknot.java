package programs;

import meshi.energy.EnergyCreator;
import meshi.energy.EvaluationException;
import meshi.energy.pairwiseNonBondedTerms.ExcludedVol.ExcludedVolCreator;
import meshi.energy.simpleEnergyTerms.angle.AngleCreator;
import meshi.energy.simpleEnergyTerms.bond.BondCreator;
import meshi.energy.simpleEnergyTerms.outOfPlane.OutOfPlaneCreator;
import meshi.energy.simpleEnergyTerms.plane.PlaneCreator;
import meshi.energy.stretch.StretchCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.sequences.AlignmentException;
import meshi.util.*;
import meshi.util.file.MeshiWriter;

import java.io.IOException;

/**
 * Created with IntelliJ IDEA.
 * User: chen
 * Date: 30/04/14
 * Time: 17:27
 * To change this template use File | Settings | File Templates.
 */
public class Unknot extends MeshiProgram{
    private static EnergyCreator[] energyCreators1 = {
            new BondCreator(),
            new AngleCreator(),
            new PlaneCreator(),
            new OutOfPlaneCreator(),
            new ExcludedVolCreator()
    };
    private static EnergyCreator[] energyCreators2 = {
            new StretchCreator(),
            new BondCreator("simple"),
            new AngleCreator(),
            new PlaneCreator(),
            new OutOfPlaneCreator(),
            new ExcludedVolCreator()
    };

    public static void main(String[] args) throws EvaluationException, UpdateableException, AlignmentException,IOException{
        initRandom(0);
        CommandList commands = new CommandList(args, new CommandsException("This is weird"));
        Protein  model = Utils.getProtein(commands,args[1],new ResidueExtendedAtomsCreator(),Utils.defaultExceptionHandler);
        model.atoms().sideChains().setNowhere();
        Utils.relax(model.atoms(), model, energyCreators1, commands);
        Utils.relax(model.atoms(),model,energyCreators2, commands);
        MeshiWriter writer = new MeshiWriter(args[2]);
        model.atoms().located().print(writer);
    }
}
