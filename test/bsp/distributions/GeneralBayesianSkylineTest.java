package bsp.distributions;


import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.BayesianSkyline;
import beast.evolution.tree.coalescent.TreeIntervals;
import beast.util.TreeParser;
import junit.framework.TestCase;
import org.junit.Test;


/**
 * @author Louis du Plessis
 * @date   2018/09/16
 */
public class GeneralBayesianSkylineTest extends TestCase {


    @Test
    public void test1 () {

        Tree tree = new TreeParser("((D4Philip56:2.0,(D4Philip64:3.0,D4Philip84:23.0):7.0):10.0,(D4SLanka78:19.0,(D4Thai78:5.0,D4Thai84:11.0):14.0):15.0);",false);
        TreeIntervals intervals = new TreeIntervals(tree);

        GeneralBayesianSkyline skyline = new GeneralBayesianSkyline();
        skyline.initByName("popSizes", "1.0 2.0", "groupSizes", "2 3", "treeIntervals", intervals);

        BayesianSkyline skyline2 = new BayesianSkyline();
        skyline2.initByName("popSizes", "1.0 2.0", "groupSizes", "2 3", "treeIntervals", intervals);

        double logP = skyline.calculateLogP();
        double logP2 = skyline2.calculateLogP();

        System.out.println(logP+"\t"+logP2);

    }

    @Test
    public void test2 () {

        Tree tree = new Tree("(((1:1,2:1):2.5,(3:1.5,4:1.5):2):2,5:5.5);");
        TreeIntervals intervals = new TreeIntervals(tree);

        GeneralBayesianSkyline skyline = new GeneralBayesianSkyline();
        skyline.initByName("popSizes", "1.0 2.0", "groupSizes", "2 2", "treeIntervals", intervals);

        BayesianSkyline skyline2 = new BayesianSkyline();
        skyline2.initByName("popSizes", "1.0 2.0", "groupSizes", "2 2", "treeIntervals", intervals);

        double logP = skyline.calculateLogP();
        double logP2 = skyline2.calculateLogP();

        System.out.println(logP+"\t"+logP2);

    }

    public void test3 () {

        Tree tree = new TreeParser("((D4Philip56:2.0,(D4Philip64:3.0,D4Philip84:23.0):7.0):10.0,(D4SLanka78:17.0,(D4Thai78:5.0,D4Thai84:11.0):12.0):17.0);",false);
        TreeIntervals intervals = new TreeIntervals(tree);

        GeneralBayesianSkyline skyline = new GeneralBayesianSkyline();
        skyline.initByName("popSizes", "1.0 2.0", "groupSizes", "2 3", "treeIntervals", intervals);

        BayesianSkyline skyline2 = new BayesianSkyline();
        skyline2.initByName("popSizes", "1.0 2.0", "groupSizes", "2 3", "treeIntervals", intervals);

        double logP = skyline.calculateLogP();
        double logP2 = skyline2.calculateLogP();

        System.out.println(logP+"\t"+logP2);

    }


}
