package bsp.distributions;

import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.TreeIntervals;
import beast.util.TreeParser;
import org.junit.Test;

public class PreferentialBayesianSkylineTest {

    @Test
    public void test1 () {

        Tree tree = new TreeParser("((D4Philip56:2.0,(D4Philip64:3.0,D4Philip84:23.0):7.0):10.0,(D4SLanka78:19.0,(D4Thai78:5.0,D4Thai84:11.0):14.0):15.0);",false);
        TreeIntervals intervals = new TreeIntervals(tree);

        PreferentialBayesianSkyline skyline = new PreferentialBayesianSkyline();
        skyline.initByName("popSizes", "1.0 2.0 3.0", "groupSizes", "4 4 3", "treeIntervals", intervals, "samplingIntensity", "0.4");


        double logP = skyline.calculateLogP();

        System.out.println(logP);

    }
}
