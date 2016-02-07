package meshi.scoringFunctions;

import meshi.energy.EvaluationException;
import meshi.energy.TotalEnergy;
import meshi.molecularElements.Protein;
import meshi.util.Command;
import meshi.util.CommandList;
import meshi.util.UpdateableException;
import meshi.util.info.DoubleInfoElementList;
import org.xml.sax.SAXException;

import javax.xml.parsers.ParserConfigurationException;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: chen
 * Date: 02/02/13
 * Time: 14:47
 * To change this template use File | Settings | File Templates.
 */
public class GdtScore extends CombinedEnergyScore  {
    public GdtScore(String parametersFileName, String name)throws IOException,EvaluationException, UpdateableException,ParserConfigurationException,SAXException{
        super(parametersFileName, name);
        }
    public double score(DoubleInfoElementList energyInfo) throws EvaluationException, UpdateableException{
        double temp = super.score(energyInfo);
        if (temp >= 1) return 1;
        if (temp <= 0.11) return 0.11;
        return temp;
    }

    public static ArrayList<Score> getScoreFunctions(CommandList commands) throws IOException, UpdateableException, EvaluationException,ParserConfigurationException, SAXException{
        CommandList scoreCommands = commands.firstWordFilter("selectionScore");
        if (scoreCommands.size() == 0) return null;
        ArrayList<Score> out = new ArrayList<Score>(scoreCommands.size());
        for (Command command:scoreCommands) {
            out.add(new GdtScore(command.secondWord(),command.thirdWord()));
        }
        return out;
    }
}
