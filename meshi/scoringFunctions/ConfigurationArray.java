package meshi.scoringFunctions;

import org.w3c.dom.Document;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Created with IntelliJ IDEA.
 * User: chen
 * Date: 11/12/13
 * Time: 19:19
 * To change this template use File | Settings | File Templates.
 */
public class ConfigurationArray extends ArrayList<Configuration>{
        public ConfigurationArray(File xmlFile) throws IOException,ParserConfigurationException , SAXException{
            DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
            DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
            Document doc = dBuilder.parse(xmlFile);
            doc.getDocumentElement().normalize();
            NodeList nodes = doc.getElementsByTagName("configuration");
            System.out.println("xxxxxxxxxx "+nodes.getLength());
            if (nodes.getLength() == 0)
                throw new RuntimeException("failed to read file "+xmlFile);
            for (int iNode = 0; iNode < nodes.getLength(); iNode++) {
                add(new Configuration(nodes.item(iNode)));
            }
        }

}
