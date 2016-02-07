/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package programs;

import meshi.energy.EnergyCreator;
import meshi.energy.EnergyInfoElement;
import meshi.energy.EvaluationException;
import meshi.energy.TotalEnergy;
import meshi.energy.beta.BetaCreator;
import meshi.energy.conservationContacts.ConservationContactsCreator;
import meshi.energy.conservationContacts.ConservationContactsHrCreator;
import meshi.energy.hydrogenBond.HydrogenBondsCreator;
import meshi.energy.hydrogenBondsAngle.HBondsPunishOHNAngleCreator;
import meshi.energy.hydrogenBondsAngle.HbondsPunishHOCAngleCreator;
import meshi.energy.hydrogenBondsPairs.HydrogenBondsPairsCreator;
import meshi.energy.one.OneCreator;
import meshi.energy.pairwiseNonBondedTerms.CooperativeSumma.CooperativeZ5typesSummaCreator;
import meshi.energy.pairwiseNonBondedTerms.CooperativeSumma.CooperativeZStd5typesSummaCreator;
import meshi.energy.pairwiseNonBondedTerms.ExcludedVol.ExcludedVolCreator;
import meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma.AtomicPairwisePMFSumma;
import meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma.AtomicPairwisePMFSummaCreator;
import meshi.energy.rapdf.RapdfCreator;
import meshi.energy.rg.RgCreator;
import meshi.energy.secondaryStructure.SecondaryStructureFractionCreator;
import meshi.energy.simpleEnergyTerms.DistanceConstraints.DistanceConstraintsCreator;
import meshi.energy.simpleEnergyTerms.angle.AngleCreator;
import meshi.energy.simpleEnergyTerms.bond.BondCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ResidueTorsions;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ResidueTorsionsList;
import meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity.CompositePropensityCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity.CooperativeZPropensityCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity.CooperativeZStdPropensityCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandran.RamachandranCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandran.RamachandranEnergy;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandran.RamachandranEnergyElement;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranCore.RamachandranCoreCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.CooperativeZRamachandranCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.CooperativeZStdRamachandranCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergyCreator;
import meshi.energy.simpleEnergyTerms.inflate.inflate.InflateCreator;
import meshi.energy.simpleEnergyTerms.inflate.inflate.InflateType;
import meshi.energy.simpleEnergyTerms.outOfPlane.OutOfPlaneCreator;
import meshi.energy.simpleEnergyTerms.plane.PlaneCreator;
import meshi.energy.simpleEnergyTerms.plane.PlaneEnergy;
import meshi.energy.simpleEnergyTerms.plane.PlaneParametersList;
import meshi.energy.simpleEnergyTerms.tether.TetherCreator;
import meshi.energy.simpleEnergyTerms.tether.TetherEnergy;
import meshi.energy.solvate.SolvateCreatorHBforMinimization;
import meshi.energy.twoTorsions.FlatRamachCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.ConSeq;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.ResidueList;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtoms;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.optimizers.*;
import meshi.parameters.MeshiPotential;
import meshi.parameters.SecondaryStructure;
import meshi.scoringFunctions.CombinedEnergyScore;
import meshi.scoringFunctions.GdtScore;
import meshi.scoringFunctions.EnergyScore;
import meshi.scoringFunctions.Score;
import meshi.sequences.AlignmentException;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.util.*;
import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;
import meshi.util.filters.Filter;
import meshi.util.info.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;

public class Optimize extends MeshiProgram implements KeyWords {
    private enum TETHER_FLAG {
        ON, OFF, RESET_ON, RESET_OFF
    }

    private enum SUMMA_COOP_FLAG {
        ON, OFF
    }



    private enum RAMACH_COOP_FLAG {
        ON, OFF
    }



    private enum PROP_COOP_FLAG {
        ON, OFF
    }



    private enum RG_FLAG {
        ON, OFF
    }

    private enum SOLVATE_FLAG{
         ON, OFF
    }


    public static final String NAME = "Optimize";

    private static MCM mcm;
    private static Relaxation relaxation;
    private static LBFGS lbfgs;
    private static String inFileName, nativeFileName, outFileName,dsspFile;
    private static Protein model, originalModel;
    private static Boolean OK;
    private static CommandList commands;
    private static TotalEnergy minimizationEnergy;
    private static ArrayList<Score> scoreFunctions;
    private static TotalEnergy perturbationEnergy1, perturbationEnergy2, perturbationEnergy3, perturbationEnergy4;
    private static DistanceMatrix distanceMatrix;
    private static OptimizeLogger log;
    private static int seed;


    private static TetherCreator tetherAllCreator = new TetherCreator(InfoType.TETHER_ALL);
    private static TetherCreator residueTetherCreator = new TetherCreator(InfoType.TETHER);
    private static HydrogenBondsCreator hydrogenBondsCreator = new HydrogenBondsCreator();
    private static SolvateCreatorHBforMinimization solvateCreator = new SolvateCreatorHBforMinimization();
    private static RgCreator rgCreator = new RgCreator();
    private static RapdfCreator rapdfCreator = new RapdfCreator();
    private static RapdfCreator rapdfCreator2 = new RapdfCreator("Rapdf2", "rapdfoutputfile2");

    private static CompositePropensityCreator propensityCreator = new CompositePropensityCreator();
    private static CooperativeZPropensityCreator cooperativeZPropensityCreator = new CooperativeZPropensityCreator(propensityCreator);
    private static CooperativeZStdPropensityCreator cooperativeZStdPropensityCreator = new CooperativeZStdPropensityCreator(cooperativeZPropensityCreator);

    private static AtomicPairwisePMFSummaCreator summaCreator = new AtomicPairwisePMFSummaCreator();
    private static CooperativeZ5typesSummaCreator cooperativeZSummaCreator = new CooperativeZ5typesSummaCreator(summaCreator);
    private static CooperativeZStd5typesSummaCreator cooperativeZStdSummaCreator = new CooperativeZStd5typesSummaCreator(cooperativeZSummaCreator);

    private static ExcludedVolCreator excludedVolCreator = new ExcludedVolCreator();
    private static RamachandranSidechainEnergyCreator ramachCreator = new RamachandranSidechainEnergyCreator();
    private static RamachandranCoreCreator ramachandranCoreCreator = new RamachandranCoreCreator(ramachCreator);
    private static CooperativeZRamachandranCreator cooperativeZRamachandranCreator = new CooperativeZRamachandranCreator(ramachCreator);
    private static CooperativeZStdRamachandranCreator cooperativeZStdRamachandranCreator = new CooperativeZStdRamachandranCreator(cooperativeZRamachandranCreator);
    private static RamachandranCreator ramachandranCreator = new RamachandranCreator();
    private static BetaCreator betaCreator = new BetaCreator();
    private static InflateCreator inflateCreator = new InflateCreator(InflateType.SIMPLE);
    private static InflateCreator inflatePerSegmentCreator = new InflateCreator(InflateType.PER_SEGMENT);
    private static InflateCreator inflateBySegmentCreator = new InflateCreator(InflateType.BY_SEGMENT);
    private static InflateCreator inflateByOtherModelCreator = new InflateCreator(InflateType.BY_OTHER_MODEL, new File("."), new MyFileFilter());
    private static ConservationContactsCreator conservationContactsCreator8 = new ConservationContactsCreator(InfoType.CONSERVATION_CONTACTS8);
    private static ConservationContactsCreator conservationContactsCreator11 = new ConservationContactsCreator(InfoType.CONSERVATION_CONTACTS11);
    private static ConservationContactsCreator conservationContactsCreator15 = new ConservationContactsCreator(InfoType.CONSERVATION_CONTACTS15);
    private static ConservationContactsHrCreator conservationContactsHrCreator = new ConservationContactsHrCreator();
    private static HydrogenBondsPairsCreator hydrogenBondsPairsCreator = new HydrogenBondsPairsCreator(hydrogenBondsCreator);
    private static  FlatRamachCreator flatRamachCreator = new FlatRamachCreator();    //private static ExcludedVolCreator excludedVolCreator = new ExcludedVolCreator();
    private static OneCreator oneCreator = new OneCreator();
    private static PlaneCreator planeCreator = new PlaneCreator();
    private Score optimizationScore = null;
    private static EnergyCreator[] energyCreators = {
            new DistanceConstraintsCreator(),
            new BondCreator(),
            new AngleCreator(),
            planeCreator,
            new OutOfPlaneCreator(),
            //excludedVolCreator,
            ramachandranCreator,
            ramachCreator,
            ramachandranCoreCreator,
            cooperativeZRamachandranCreator,
            cooperativeZStdRamachandranCreator,
            propensityCreator,
            cooperativeZPropensityCreator,
            cooperativeZStdPropensityCreator,
            summaCreator,
            excludedVolCreator,
            cooperativeZSummaCreator,
            cooperativeZStdSummaCreator,
            solvateCreator,
            hydrogenBondsCreator,
            hydrogenBondsPairsCreator,
            new HbondsPunishHOCAngleCreator(hydrogenBondsCreator),
            new HBondsPunishOHNAngleCreator(hydrogenBondsCreator),
            rgCreator,
            conservationContactsCreator8,
            conservationContactsCreator11,
            conservationContactsCreator15,
            conservationContactsHrCreator,
            oneCreator,
            betaCreator,
            inflateCreator,
            inflatePerSegmentCreator,
            inflateBySegmentCreator,
            inflateByOtherModelCreator,
            residueTetherCreator,
            flatRamachCreator,
                tetherAllCreator,
						     rapdfCreator,
              rapdfCreator2,
            new SecondaryStructureFractionCreator()

    };
    private static EnergyCreator[] excludeFromMinimization = {betaCreator, inflateCreator, inflatePerSegmentCreator,
            inflateBySegmentCreator, inflateByOtherModelCreator, ramachandranCreator};       // inflate is noise  and the other two have almost arbitrary values.
    private static EnergyCreator[] excludeFromPerturbation1 = {tetherAllCreator, inflatePerSegmentCreator, inflateBySegmentCreator,
                                                               inflateByOtherModelCreator, rapdfCreator, rapdfCreator2};
    private static EnergyCreator[] excludeFromPerturbation2 = {tetherAllCreator, inflateCreator, inflateBySegmentCreator,
                                                               inflateByOtherModelCreator , rapdfCreator,  rapdfCreator2};
    private static EnergyCreator[] excludeFromPerturbation3 = {tetherAllCreator, inflateCreator, inflatePerSegmentCreator,
                                                               inflateByOtherModelCreator , rapdfCreator, rapdfCreator2};
    private static EnergyCreator[] excludeFromPerturbation4 = {tetherAllCreator, inflateCreator, inflatePerSegmentCreator,
                                                               inflateBySegmentCreator, rapdfCreator, rapdfCreator2};
    private static String parentString;

    public static void main(String[] argv) throws IOException, OptimizerException, UpdateableException, EvaluationException {
        init(argv);
        // Chemistry
        model = null;
        try {
            parentString = getParentString(inFileName);
            model = Utils.getProtein(commands, inFileName, ResidueExtendedAtomsCreator.creator, Utils.defaultExceptionHandler);
            conservationContactsCreator8.setFileFound(ConSeq.setConSeq(model.chain(), inFileName));
            conservationContactsCreator11.setFileFound(ConSeq.setConSeq(model.chain(), inFileName));
            conservationContactsCreator15.setFileFound(ConSeq.setConSeq(model.chain(), inFileName));
            conservationContactsHrCreator.setFileFound(ConSeq.setConSeq(model.chain(), inFileName));
            originalModel = Utils.getProtein(commands, inFileName, ResidueExtendedAtomsCreator.creator, Utils.defaultExceptionHandler);
            OK = Utils.addAtoms(model, false, commands,
                                new PlaneParametersList(EnergyCreator.parametersDirectory(commands) +
                                                                           "/" + MeshiPotential.PLANE_PARAMETERS),false);
            Utils.println("nAtoms " + model.atoms().size());
            Utils.println("nResidues " + model.residues().size());
            if (!OK) throw new RuntimeException("Too many atoms are missing. Cannot refine");
            ((MyFileFilter) inflateByOtherModelCreator.filter).reset(model);
            model.atoms().defrost();

    //      Secondary structure;
            if (!dsspFile.equals("NONE"))
                Utils.AssignDSSP(model, dsspFile);
            Utils.setSS(model, commands);
            Utils.println("secondary structure 1: ");
            int i = 0;
            for (Residue residue : model.residues()) {
                Utils.print(residue.secondaryStructure() + "  ");
                if (i++ % 10 == 0) Utils.println();
            }
            Utils.println();

            // conserved residues
            ResidueList conservedResidues = getConservedResidues(model, commands);

            // Work power
            perturbationEnergy1 = new TotalEnergy(model,
                    excludeEnergy(energyCreators, excludeFromPerturbation1),
                    commands,"perturbationEnergy1");
            perturbationEnergy2 = new TotalEnergy(model,
                    excludeEnergy(energyCreators, excludeFromPerturbation2),
                    commands,"perturbationEnergy2");
            perturbationEnergy3 = new TotalEnergy(model,
                    excludeEnergy(energyCreators, excludeFromPerturbation3),
                    commands,"perturbationEnergy3");
            perturbationEnergy4 = new TotalEnergy(model,
                    excludeEnergy(energyCreators, excludeFromPerturbation4),
                    commands,"perturbationEnergy4");
            minimizationEnergy = new TotalEnergy(model,
                    excludeEnergy(energyCreators, excludeFromMinimization),
                    commands,"minimization energy");
            ((TetherEnergy) residueTetherCreator.term()).setResidue(commands);
            scoreFunctions = GdtScore.getScoreFunctions(commands);
            ArrayList<Score> energyScores= new ArrayList<Score>();
            for (Score scoreFunction : scoreFunctions)
                if (scoreFunctions.toString().equals("energy"))
                    energyScores.add(scoreFunction);
            minimizationEnergy.setCurrent();
            //((PlaneEnergy)planeCreator.term()).fixCisTrans();

            lbfgs = Utils.getLBFGS(minimizationEnergy, commands, RELAX);
            log = new OptimizeLogger(model, originalModel, nativeFileName, outFileName, minimizationEnergy, parentString);

            //log.log("BEGINNING", true);
            log.mcm(scoreFunctions, minimizationEnergy, "BEGINNING");
            if (commands.firstWordFilter("MCM").secondWord("maxSteps").thirdWord().equals("0")){
                System.out.println(" Simulation of 0 steps done.");
                return;
            }
            minimizeCassette(lbfgs, minimizationEnergy, (TetherEnergy) tetherAllCreator.term());

            log.setEnergy(minimizationEnergy);
//            relaxation = getRelaxation(model, minimizationEnergy, scoreForRelaxation, selectionScore, perturbationEnergy1, perturbationEnergy2, perturbationEnergy3, perturbationEnergy3, commands, conservedResidues);
//            Optimizer.OptimizerStatus status = relaxation.run(log);
//            if ((status == Optimizer.OptimizerStatus.DONE) | (status == Optimizer.OptimizerStatus.CONVERGED)) {
                mcm = getMCM(model, minimizationEnergy, scoreFunctions, perturbationEnergy1, perturbationEnergy2, perturbationEnergy3, perturbationEnergy3, commands, conservedResidues);
                mcm.run(log);
//            }

        tetherAllCreator.term().off();
        log.mcm(scoreFunctions, minimizationEnergy, "MCM_END\" step=\""+mcm.maxSteps);
            mcm.energy().distanceMatrix().writeNonBondedList(mcm.maxSteps+"_MCM_END");

        }
         catch (Exception ex){
             MeshiWriter writer = new MeshiWriter(outFileName);
             writer.println("Simulation failed due to:\n"+ex+"\n");
             ex.printStackTrace();
             if (model  == null)
                writer.println("No model was created.");
             else if (model.atoms().size() <= 0)
                 writer.println("No atoms were created in the model.");
                 else model.atoms().print(writer);
             writer.close();
             throw new RuntimeException("Simulation failed.");
         }

    }


/*    private static TotalEnergy perturbationEnergy1() throws UpdateableException{
        return new TotalEnergy(model,
                excludeEnergy(energyCreators, excludeFromPerturbation1),
                commands,"perturbationEnergy1");
    }
    private static TotalEnergy perturbationEnergy2()throws UpdateableException{
        return new TotalEnergy(model,
            excludeEnergy(energyCreators, excludeFromPerturbation2),
            commands,"perturbationEnergy2");
    }
    private static TotalEnergy perturbationEnergy3() throws UpdateableException {
        return new TotalEnergy(model,
                excludeEnergy(energyCreators, excludeFromPerturbation3),
                commands,"perturbationEnergy3");
    }
    private static TotalEnergy perturbationEnergy4() throws UpdateableException{
        return new TotalEnergy(model,
            excludeEnergy(energyCreators, excludeFromPerturbation4),
            commands,"perturbationEnergy4");
    }
    private static TotalEnergy minimizationEnergy() throws UpdateableException{
        return new TotalEnergy(model,
            excludeEnergy(energyCreators, excludeFromMinimization),
            commands,"minimization energy");
    }
*/
    private static ResidueList getConservedResidues(Protein model, CommandList commands) {
        ResidueList out = new ResidueList();
        String line = commands.firstWord("conserved").secondWord();
        String[] words = line.split(",");
        for (String word : words) {
            int number = Integer.valueOf(word);
            for (Residue residue : model.residues()) {
                if (residue.number() == number) {
                    out.add(residue);
                    Utils.println(residue + " marked as conserved");
                }
            }
        }
        return out;
    }

    private static void refreshBeta(RamachandranEnergy ramachandranEnergy) {
        ResidueTorsions rt1 = null, rt2 = null;
        ResidueTorsionsList residueTorsionsList = ramachandranEnergy.residueTorsionsList();
        for (ResidueTorsions residueTorsions : residueTorsionsList) {
            rt2 = residueTorsions;
            if (rt1 != null) {
                if (rt1.isBeta() && rt2.isBeta()) {
                    Residue residue = rt2.residue();
                    SecondaryStructure secondaryStructure = residue.secondaryStructure();
                    if (!secondaryStructure.equals(SecondaryStructure.SHEET))
                        Utils.println("Changing the secondary structure of " + residue + " to SHEET.");
                    rt2.residue().setSecondaryStructure(SecondaryStructure.SHEET);
                }
            }
            rt1 = rt2;
        }
    }

    private static ArrayList<RamachandranEnergyElement> ramachEnergyPerResidue(RamachandranEnergy ramach) {
        ArrayList<RamachandranEnergyElement> out = new ArrayList<RamachandranEnergyElement>();
        for (Object o : ramach.elementsList()) {
            RamachandranEnergyElement element = (RamachandranEnergyElement) o;
            double e = element.evaluate() / element.weight();
            if (e > 5) Utils.println("Ramachandran element " + element + " ; minimizationEnergy = " + e);
            if (e > RamachandranEnergyElement.THRESHOLD) out.add(element);
        }
        return out;
    }


    private static MCM getMCM(Protein model,
                                  TotalEnergy minimizationEnergy,
                                  ArrayList<Score> scoreFunctions,
                                  TotalEnergy perturbationEnergy1,
                                  TotalEnergy perturbationEnergy2,
                                  TotalEnergy perturbationEnergy3,
                                  TotalEnergy perturbationEnergy4,
                                  CommandList commands,
                                  ResidueList conservedResidues) throws UpdateableException,EvaluationException {
        Score optimizationScore = null;
        for (Score score:scoreFunctions)
            if (score.toString().equals("optimizationScore"))  {
                optimizationScore = score;
                break;
            }
        if (optimizationScore == null)
            throw new RuntimeException("No optimization score.");

        return getMCM(model,minimizationEnergy,scoreFunctions,optimizationScore,
                                         perturbationEnergy1,perturbationEnergy2,perturbationEnergy3,perturbationEnergy4,
                                         commands, conservedResidues, MCM.mcmMode.OPTIMIZATION);
        }
 /*   private static Relaxation getRelaxation(Protein model,
                                  TotalEnergy minimizationEnergy,
                                  ArrayList<Score> scoreFunctions,
                                  TotalEnergy perturbationEnergy1,
                                  TotalEnergy perturbationEnergy2,
                                  TotalEnergy perturbationEnergy3,
                                  TotalEnergy perturbationEnergy4,
                                  CommandList commands,
                                  ResidueList conservedResidues) throws UpdateableException,EvaluationException {
            return (Relaxation) getMCM(model, minimizationEnergy, scoreFunctions,
                                                                  perturbationEnergy1, perturbationEnergy2,perturbationEnergy3,perturbationEnergy4,
                                                                  commands, conservedResidues, MCM.mcmMode.RELAXATION);
        }
   */
    private static MCM getMCM(Protein model,
                              TotalEnergy minimizationEnergy,
                              ArrayList<Score> scoreFunctions,
                              Score optimizationScore,
                              TotalEnergy perturbationEnergy1,
                              TotalEnergy perturbationEnergy2,
                              TotalEnergy perturbationEnergy3,
                              TotalEnergy perturbationEnergy4,
                              CommandList commands,
                              ResidueList conservedResidues, MCM.mcmMode mode) throws UpdateableException,EvaluationException {
        Perturbation perturbation1 = new PerturbationByMinimization(perturbationEnergy1, commands, model, conservedResidues, " perturbation with a complete energy function");
        Perturbation perturbation2 = new PerturbationByMinimization(perturbationEnergy2, commands, model, conservedResidues, " perturbation with a minimal energy function");
        Perturbation perturbation3 = new PerturbationByMinimization(perturbationEnergy3, commands, model, conservedResidues, " perturbation with a minimal energy function");
        Perturbation perturbation4 = new PerturbationByMinimization(perturbationEnergy4, commands, model, conservedResidues, " perturbation with a minimal energy function");
        Perturbation[] perturbations = {perturbation1, perturbation2, perturbation3};
        Perturbation perturbation = new CombinedPerturbation(perturbations);
        //new ScmodPerturbations(model, commands,conservedResidues)};
        minimizationEnergy.setCurrent();
        Minimizer minimizer = Utils.getLBFGS(minimizationEnergy, commands, MINIMIZE);
        double initialTemperature = commands.firstWordFilter(MC_MINIMIZATION).secondWord(INITIAL_TEMPERATURE).thirdWordDouble();
        double finalTemperature = commands.firstWordFilter(MC_MINIMIZATION).secondWord(FINAL_TEMPERATURE).thirdWordDouble();
        int nSteps = commands.firstWordFilter(MC_MINIMIZATION).secondWord(MAX_STEPS).thirdWordInt();
        TemperatureGenerator temperatureGenerator = new TemperatureGenerator(initialTemperature, finalTemperature, nSteps);
        //AbstractEnergy[] excludedTerms = {inflateCreator.term(), tetherAllCreator.term()};
        if (mode == MCM.mcmMode.RELAXATION )
            return new Relaxation(minimizationEnergy, scoreFunctions, optimizationScore, minimizer, perturbation, temperatureGenerator, nSteps);
        return new MCM(minimizationEnergy, scoreFunctions, optimizationScore, minimizer, perturbation, temperatureGenerator, nSteps);
    }

    private static Optimizer.OptimizerStatus minimize(LBFGS lbfgs, TotalEnergy energy,
                                 TETHER_FLAG tetherFlag,
                                 RG_FLAG rgFlag,SOLVATE_FLAG solvateFlag) throws OptimizerException, UpdateableException,EvaluationException {

        return minimize(lbfgs, energy, RAMACH_COOP_FLAG.ON,
                PROP_COOP_FLAG.ON,
                SUMMA_COOP_FLAG.ON,
                tetherFlag,
                rgFlag,solvateFlag);
    }

    private static Optimizer.OptimizerStatus minimize(LBFGS lbfgs,TotalEnergy energy,
                                 RAMACH_COOP_FLAG rcFlag,
                                 PROP_COOP_FLAG prFlag,
                                 SUMMA_COOP_FLAG scFlag,
                                 TETHER_FLAG tetherFlag,
                                 RG_FLAG rgFlag,
                                 SOLVATE_FLAG solvateFlag) throws OptimizerException, UpdateableException,EvaluationException {


        if ((rcFlag == RAMACH_COOP_FLAG.OFF) &&
            (cooperativeZRamachandranCreator.term() != null)) {
            cooperativeZRamachandranCreator.term().off();
            cooperativeZStdRamachandranCreator.term().off();
        }
        if ((prFlag == PROP_COOP_FLAG.OFF) &&
            (cooperativeZPropensityCreator.term() != null)){
            cooperativeZPropensityCreator.term().off();
            cooperativeZStdPropensityCreator.term().off();
        }
        if ((scFlag == SUMMA_COOP_FLAG.OFF) &&
            (cooperativeZSummaCreator.term()!= null)) {
            cooperativeZSummaCreator.term().off();
            cooperativeZStdSummaCreator.term().off();
        }
        if (tetherFlag == TETHER_FLAG.OFF) {
            tetherAllCreator.term().off();
        } else if (tetherFlag == TETHER_FLAG.RESET_OFF) {
            ((TetherEnergy) tetherAllCreator.term()).reset();
            tetherAllCreator.term().off();
        } else if (tetherFlag == TETHER_FLAG.RESET_ON) {
            ((TetherEnergy) tetherAllCreator.term()).reset();
        }
        if  (rgCreator.term() != null) {
            if(rgFlag == RG_FLAG.OFF) rgCreator.term().off();
            else rgCreator.term().on();
        }
        if (solvateFlag == SOLVATE_FLAG.OFF) {
            solvateCreator.term().off();
        }
        RamachandranEnergy.resetRamachandran(energy);

        return lbfgs.run();
    }

    private static void minimizeCassette(LBFGS lbfgs, TotalEnergy energy, TetherEnergy tether) throws OptimizerException, IOException, UpdateableException, EvaluationException,AlignmentException {
        Optimizer.OptimizerStatus os = null;
        Utils.println("\nStarting first phase of minimizationCassette \n");

                rgCreator.term().scaleWeight(0.1);
                tether.scaleWeight(0.1);
                ramachCreator.term().scaleWeight(0.1);
                hydrogenBondsPairsCreator.term().scaleWeight(0.1);
                ((PlaneEnergy)planeCreator.term()).scaleWeight(0.01);
                Utils.print("minimizationCassete   1");
                minimize(lbfgs,energy,
                        RAMACH_COOP_FLAG.OFF,
                        PROP_COOP_FLAG.OFF,
                        SUMMA_COOP_FLAG.OFF,
                        TETHER_FLAG.RESET_ON,
                        RG_FLAG.ON,SOLVATE_FLAG.ON);

                log.log("REFINE 0 ", Utils.verbose());
                rgCreator.term().scaleWeight(2);
                tether.scaleWeight(10);
                ramachCreator.term().scaleWeight(10);
                hydrogenBondsPairsCreator.term().scaleWeight(10);
                ((PlaneEnergy)planeCreator.term()).scaleWeight(10);
                Utils.print("minimizationCassete   2");
                minimize(lbfgs,energy,
                        RAMACH_COOP_FLAG.OFF,
                        PROP_COOP_FLAG.OFF,
                        SUMMA_COOP_FLAG.OFF,
                        TETHER_FLAG.RESET_ON,
                        RG_FLAG.ON,SOLVATE_FLAG.ON);
                log.log("REFINE 0.1 ", Utils.verbose());
                rgCreator.term().scaleWeight(5);
                ((PlaneEnergy)planeCreator.term()).scaleWeight(10);
                    Utils.print("minimizationCassete   3.");
                    minimize(lbfgs,energy,
                            RAMACH_COOP_FLAG.OFF, PROP_COOP_FLAG.OFF,
                            SUMMA_COOP_FLAG.OFF, TETHER_FLAG.RESET_ON,
                            RG_FLAG.ON, SOLVATE_FLAG.ON);
                log.log("REFINE 0.2 ", Utils.verbose());

        if (cooperativeZSummaCreator.term() != null)
            cooperativeZSummaCreator.term().on();
        double eTest = energy.evaluate();
        Utils.println("Energy with cooperativeSumma = "+eTest);
    //   while (((AtomicPairwisePMFSumma) summaCreator.term()).numberOfClashes()>0) {
        while (eTest > 100000) {
            int nClashes = ((AtomicPairwisePMFSumma) summaCreator.term()).numberOfClashes();
            Utils.println("\n"+nClashes+"clashes\n");
            if (nClashes > 0) Utils.println(((AtomicPairwisePMFSumma) summaCreator.term()).clashes());

                excludedVolCreator.term().scaleWeight(2);

                minimize(lbfgs,energy,
                        RAMACH_COOP_FLAG.ON, PROP_COOP_FLAG.ON,
                        SUMMA_COOP_FLAG.OFF, TETHER_FLAG.RESET_ON,
                        RG_FLAG.ON,SOLVATE_FLAG.ON);
            if (cooperativeZSummaCreator.term() != null)
                cooperativeZSummaCreator.term().on();
            eTest = energy.evaluate();
            Utils.println("Energy with cooperativeSummaPotential = "+eTest);
        }
        if (excludedVolCreator.term() != null)
            excludedVolCreator.term().off();
                log.log("REFINE 2 ", Utils.verbose());
        Utils.println("\nStarting second phase of minimizationCassette \n");
        try {
            for (int i = 0; (i < 10) & (os != Optimizer.OptimizerStatus.CONVERGED) ; i++) {
                System.out.println(energy.reportHeader());
                System.out.println(energy.report(22200));
                Utils.print("REFINE 2."+i) ;
                os = minimize(lbfgs, energy,
                    TETHER_FLAG.RESET_ON,
                    RG_FLAG.ON,SOLVATE_FLAG.ON);
                 log.log("REFINE 2."+i, Utils.verbose());
            }
            os = minimize(lbfgs, energy,
                    TETHER_FLAG.OFF, RG_FLAG.ON,SOLVATE_FLAG.ON);
            log.log("REFINE 3 ", Utils.verbose());

        }
        catch (OptimizerException ex) {
            Utils.print("\nException caught\n" + ex + "\n");
            try {
                if (Utils.verbose()) energy.test();
            } catch (UpdateableException ae) {
                System.out.println("energy.test failed due to " + ae);
                ae.printStackTrace();
                throw new RuntimeException(ae+"\n"+ae.getStackTrace()+"\nquiting");
            }
            minimize(lbfgs,energy,
                    RAMACH_COOP_FLAG.OFF, PROP_COOP_FLAG.OFF,
                    SUMMA_COOP_FLAG.ON, TETHER_FLAG.RESET_ON,
                    RG_FLAG.OFF,SOLVATE_FLAG.ON);
        }
    }

    private static class OptimizeLogger extends ProteinInfoList implements Logger {
        String outFileName, infoFileName;
        MeshiWriter output;
        MeshiInfoXMLwriter info;
        ModelAnalyzer analyzer;
        String parentString;
        TotalEnergy energy;

        public OptimizeLogger(Protein model, Protein originalModel, String nativeFileName, String outFileName, TotalEnergy energy, String parentString) {
            super("Optimization history of " + model);
            Protein nativeStructure;
            if (!(nativeFileName.equals("NONE") ||nativeFileName.equals("none"))) {
                nativeStructure = Utils.getProtein(commands, nativeFileName, ResidueExtendedAtoms.creator, Utils.defaultExceptionHandler);
            } else {
                nativeStructure = null;
            }
            analyzer = new ModelAnalyzer(model, nativeStructure, originalModel, energy, ResidueAlignmentMethod.BY_RESIDUE_NUMBER);
            this.outFileName = outFileName;
            this.parentString = parentString;
            infoFileName = outFileName + ".xml";
            this.energy = energy;
        }

        public void setEnergy(TotalEnergy energy) {
            this.energy = energy;
            analyzer.setEnergy(energy);
        }

        public double rms() {
            if(nativeFileName.equals("NONE") ||nativeFileName.equals("none") ) return -1;
            else {
                try{
                    return analyzer.rms();
                }catch (Exception ex ) {throw new RuntimeException(ex.getMessage());}
            }
        }

        public void log(String comment, boolean printFlag) throws IOException,EvaluationException,AlignmentException {

            try {
                add(analyzer.analyze(comment));
            } catch (UpdateableException ae) {
                System.out.println("log failed due to " + ae);
                ae.printStackTrace();
                throw new RuntimeException("quiting");
            }
            if (printFlag) {
                try {
                    output = new MeshiWriter(outFileName);
                    info = new MeshiInfoXMLwriter(infoFileName);
                    output.println(parentString);
                } catch (Exception ex) {
                    throw new RuntimeException(ex);
                }
                print(output);
                print(info);
                energy.on();
                try {
                    energy.evaluateAtoms();
                } catch (UpdateableException ae) {
                    System.out.println("log failed due to " + ae);
                    ae.printStackTrace();
                    throw new RuntimeException("quiting");
                } catch (EvaluationException ee) {
                    System.out.println("log failed due to " + ee);
                    ee.printStackTrace();
                    throw new RuntimeException("quiting");
                }
                Utils.colorByEnergy(model.atoms());
                analyzer.model.atoms().print(output);
                output.close();
            }
        }
         public void logXML(String comment,
                            boolean printFlag) throws IOException,EvaluationException,AlignmentException {
            try {
                add(analyzer.analyze(comment));
            } catch (UpdateableException ae) {
                System.out.println("log failed due to " + ae);
                ae.printStackTrace();
                throw new RuntimeException("quiting");
            }
            if (printFlag) {
                try {
                    output = new MeshiWriter(outFileName);
                    info = new MeshiInfoXMLwriter(infoFileName);
                    output.println(parentString);
                } catch (Exception ex) {
                    throw new RuntimeException(ex);
                }
                print(output);
                print(info);
                energy.on();
                try {
                    energy.evaluateAtoms();
                } catch (UpdateableException ae) {
                    System.out.println("log failed due to " + ae);
                    ae.printStackTrace();
                    throw new RuntimeException("quiting");
                }
                Utils.colorByEnergy(model.atoms());
                analyzer.model.atoms().print(output);
                output.close();
            }
        }


        public void mcm(ArrayList<Score> scoreFunctions, TotalEnergy energy, String label) throws IOException, UpdateableException, EvaluationException,AlignmentException {
            analyzer.setEnergy(energy);
             if (scoreFunctions == null) throw new RuntimeException("scoreFunction == null");
             long time = (new Date()).getTime() - energy.startTime;
            DoubleInfoElementList infoList = energy.energyInfo();
            for (Score scoreFunction : scoreFunctions) {
                Utils.println(" Calculating score "+scoreFunction);
//                ((CombinedEnergyScore) scoreFunction).setEnergy(energy);
                ((CombinedEnergyScore) scoreFunction).setLengthElement(model.chain().numberOfNonDummyResidues());
                double rmsFromOriginal;
                try {
                    rmsFromOriginal = Rms.rms(originalModel, model,ResidueAlignmentMethod.BY_RESIDUE_NUMBER);
                }
                catch (Exception ex) {throw new RuntimeException(ex);
                }
                ((CombinedEnergyScore) scoreFunction).setChangeElement(rmsFromOriginal);
                double score = scoreFunction.score(infoList);
                label = label + "\" "+scoreFunction.toString()+"=\""+score+" ";
            }
 	        log(label+"\" time=\""+time, true);
        }

        public void mcm(ArrayList<Score> scoreFunctions,
                        TotalEnergy energy,
                        int i,
                        MCM.mcmStepResult mcmStepResult) throws IOException,
                                                                UpdateableException,
                                                                EvaluationException,AlignmentException {
 //           mcm(optimizationScore, selectionScore, energy, "MCM\"  step=\""+i+"\" score=\""+score.score()+"\" result=\""+mcmStepResult+"\"  lastSuccess=\""+mcmStepResult.lastSuccess());
            mcm(scoreFunctions, energy, "MCM\"  step=\""+i+"\" result=\""+
                mcmStepResult+"\"  lastSuccess=\""+mcmStepResult.lastSuccess());
         }
    }


    private static void init(String[] argv) {
        //                 0            1            2            3                4              5
        String[] keys = {"commands", "inFileName", "dsspFile", "nativeFileName", "outFileName", "seed"};
        String[] arguments = getArguments(keys, argv);

        seed = Integer.parseInt(arguments[5]);
        System.out.println("seed " + seed);
        initRandom(seed);

        commands = new CommandList(arguments[0]);
        inFileName = arguments[1];
        dsspFile   = arguments[2];
        nativeFileName = arguments[3];
        outFileName = arguments[4];
        if (commands.keyExists("verbose")) Utils.verboseOn();
        else Utils.verboseOff();
    }

    private static String getParentString(String fileName) throws IOException {
        MeshiLineReader reader = new MeshiLineReader(fileName);
        String out;
        while ((out = reader.readLine()) != null) {
            if (out.startsWith("PARENT"))
                return out;
        }
        return "PARENT N/A";
    }


    private static EnergyCreator[] excludeEnergy(EnergyCreator[] source, EnergyCreator[] exclude) {
        EnergyCreator[] out = new EnergyCreator[source.length - exclude.length];
        int j = 0;
        for (int i = 0; i < source.length; i++) {
            if (notFound(source[i], exclude)) {
                out[j] = source[i];
                j++;
            }
        }
        return out;
    }

    private static boolean notFound(EnergyCreator energyCreator, EnergyCreator[] energyCreators) {
        for (int i = 0; i < energyCreators.length; i++) {
            if (energyCreators[i].equals(energyCreator)) return false;
        }
        return true;
    }

    private static class MyFileFilter implements Filter {
        private String prefix = null;
        Protein thisModel;

        public void reset(Protein thisModel) {
            this.thisModel = thisModel;
            int index = thisModel.name().indexOf('.');
            if (index != -1) prefix = thisModel.name().substring(0, index);
            else prefix = thisModel.name();

        }

        public boolean accept(Object obj) {
            if (prefix == null) throw new RuntimeException("prefix is " + prefix);
            File file = (File) obj;
            String path = file.getAbsolutePath();
            if (!path.endsWith("pdb")) return false;
            if (path.indexOf("out") == -1) return false;
            if (file.getName().startsWith(prefix)) {
                try {
                    double rms = Rms.rms(thisModel, Protein.getCAproteinFromApdbFile(file), ResidueAlignmentMethod.BY_RESIDUE_NUMBER);
                    if (rms < 1) return false;
                } catch (Exception ex ) {throw new RuntimeException(ex.getMessage());}
                return true;
            }
            return false;
        }
    }
}
