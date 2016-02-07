package programs;

import meshi.scoringFunctions.EnergyScore;
import meshi.util.file.MeshiLineReader;
import meshi.util.info.InfoType;

import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: keasar
 * Date: 31/01/11
 * Time: 15:05
 * To change this template use File | Settings | File Templates.
 */
public class Translate {
    public static void main(String[] argv) throws IOException{
        MeshiLineReader reader = new MeshiLineReader("out.txt");
        String line;
        char flag;
        while ((line = reader.readLine()) != null){
       //     System.out.println(line);
            String[] words = line.split("\t");
            flag = words[0].charAt(words[0].length()-1);
            words[0] = words[0].substring(0,words[0].length()-1);
            boolean found = false;
            EnergyScore.Normalization normalization;
            if (flag == '0') normalization = EnergyScore.Normalization.IGNORED;
            else if (flag == '1') normalization = EnergyScore.Normalization.NOT_NORMALIZED;
            else  normalization = EnergyScore.Normalization.NORMALIZED_BY_LENGTH;
            for (InfoType type : InfoType.values()){
                if (type.toString().equals(words[0]))         {
                    System.out.println("scoreWeight\t"+type.tag + "\t" + words[1]+"\t"+normalization);
                    found = true;
                }
            }
            if (! found) System.out.println("#not found: "+line);
        }
    }
}
