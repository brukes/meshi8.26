/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package programs;

import meshi.energy.EnergyCreator;
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
import meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity.CompositePropensityCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity.CooperativeZPropensityCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity.CooperativeZStdPropensityCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandran.RamachandranCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandran.RamachandranEnergy;
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
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.ResidueList;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtoms;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.optimizers.*;
import meshi.parameters.MeshiPotential;
import meshi.scoringFunctions.CombinedEnergyScore;
import meshi.scoringFunctions.EnergyScore;
import meshi.scoringFunctions.GdtScore;
import meshi.scoringFunctions.Score;
import meshi.sequences.AlignmentException;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.util.*;
import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;
import meshi.util.filters.Filter;
import meshi.util.info.InfoType;
import meshi.util.info.MeshiInfoElementList;
import meshi.util.info.MeshiInfoXMLwriter;
import meshi.util.info.ProteinInfoList;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;

public class ScoreTest extends MeshiProgram implements KeyWords {
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
    private static Score scoreForRelaxation, scoreForMCM, selectionScore;
    private static TotalEnergy perturbationEnergy1, perturbationEnergy2, perturbationEnergy3, perturbationEnergy4;
    private static DistanceMatrix distanceMatrix;
    private static OptimizeLogger log;
    private static int seed;
    private static MeshiInfoElementList energyScoreWeights;

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
    private static EnergyCreator[] energyCreators = {
            new DistanceConstraintsCreator(),
            new BondCreator(),
            new AngleCreator(),
            planeCreator,
            new OutOfPlaneCreator(),
            //excludedVolCreator,
            ramachandranCreator,
            ramachCreator,
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

    public static void main(String[] argv) throws Exception{
        init(argv);
        // Chemistry
        model = null;
//        try {
            parentString = getParentString(inFileName);
            model = Utils.getProtein(commands, inFileName, ResidueExtendedAtomsCreator.creator, Utils.defaultExceptionHandler);
            OK = Utils.addAtoms(model, false, 10000, commands, new PlaneParametersList(EnergyCreator.parametersDirectory(commands) +
                    "/" + MeshiPotential.PLANE_PARAMETERS),true);
            model.atoms().defrost();

    //      Secondary structure;
            if (!dsspFile.equals("NONE"))
                Utils.AssignDSSP(model, dsspFile);
            minimizationEnergy = new TotalEnergy(model,
                    excludeEnergy(energyCreators, excludeFromMinimization),
                    commands,"minimization energy");
            ((TetherEnergy) residueTetherCreator.term()).setResidue(commands);
            energyScoreWeights = EnergyScore.getWeights(commands);
            scoreForRelaxation = new EnergyScore(minimizationEnergy,  null,  model);
            scoreForMCM = new EnergyScore(minimizationEnergy, energyScoreWeights, model);
            //selectionScore = new GdtScore(minimizationEnergy, model, commands,"selectionScore");

            minimizationEnergy.setCurrent();
        log = new OptimizeLogger(model,model,"T0801-D1.BAKER-ROBETTA_TS1.out.1.pdb","out.pdb",minimizationEnergy,"XXX");
        ((CombinedEnergyScore)selectionScore).setChangeElement(log.change());
            double score = selectionScore.score(minimizationEnergy.energyInfo());
            System.out.println("score =  "+score);
        log.logXML("MCM_END",true);
//        }
//         catch (Exception ex){
//             throw new RuntimeException("Simulation failed. "+ex);
//         }

    }


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


        if (rcFlag == RAMACH_COOP_FLAG.OFF) {
            cooperativeZRamachandranCreator.term().off();
            cooperativeZStdRamachandranCreator.term().off();
        }
        if (prFlag == PROP_COOP_FLAG.OFF) {
            cooperativeZPropensityCreator.term().off();
            cooperativeZStdPropensityCreator.term().off();
        }
        if (scFlag == SUMMA_COOP_FLAG.OFF) {
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
                        RG_FLAG.ON, SOLVATE_FLAG.ON);

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
                        RG_FLAG.ON, SOLVATE_FLAG.ON);
                log.log("REFINE 0.1 ", Utils.verbose());
                rgCreator.term().scaleWeight(5);
                ((PlaneEnergy)planeCreator.term()).scaleWeight(10);
                    Utils.print("minimizationCassete   3.");
                    minimize(lbfgs,energy,
                            RAMACH_COOP_FLAG.OFF, PROP_COOP_FLAG.OFF,
                            SUMMA_COOP_FLAG.OFF, TETHER_FLAG.RESET_ON,
                            RG_FLAG.ON, SOLVATE_FLAG.ON);
                log.log("REFINE 0.2 ", Utils.verbose());

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
                        RG_FLAG.ON, SOLVATE_FLAG.ON);
            cooperativeZSummaCreator.term().on();
            eTest = energy.evaluate();
            Utils.println("Energy with cooperativeSummaPotential = "+eTest);
        }
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
                    RG_FLAG.ON, SOLVATE_FLAG.ON);
                 log.log("REFINE 2."+i, Utils.verbose());
            }
            os = minimize(lbfgs, energy,
                    TETHER_FLAG.OFF, RG_FLAG.ON, SOLVATE_FLAG.ON);
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
                    RG_FLAG.OFF, SOLVATE_FLAG.ON);
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
            if(nativeFileName.equals("NONE")) return -1;
            else {
                try{
                    return analyzer.rms();
                }catch (Exception ex ) {throw new RuntimeException(ex.getMessage());}
            }
        }
        public double change() throws Exception{
            return analyzer.change();
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


        public void mcm(ArrayList<Score> scoreFunctions,
                        TotalEnergy energy, String label) throws IOException, UpdateableException, EvaluationException,AlignmentException {
            analyzer.setEnergy(energy);
             if (selectionScore == null) throw new RuntimeException("electionScore == null");
             long time = (new Date()).getTime() - energy.startTime;
             double score = selectionScore.score(energy.energyInfo());
             double SSinterdecile = ((GdtScore)selectionScore).interdecileRange();
 	    log(label+"\" optimizationScore=\""+scoreFunctions.get(0).score(energy.energyInfo())+"\" selectionScore=\""+
                   + score+
                "\" SSinterdecile=\""+SSinterdecile+
                       "\" time=\""+time, true);
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
        String[] keys = {"commands", "inFileName", "dsspFile","seed"};
        String[] arguments = getArguments(keys, argv);

        seed = Integer.parseInt(arguments[3]);
        System.out.println("seed " + seed);
        initRandom(seed);

        commands = new CommandList(arguments[0]);
        inFileName = arguments[1];
        dsspFile   = arguments[2];
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
