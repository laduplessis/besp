package bsp.distributions;

import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.BayesianSkyline;
import beast.evolution.tree.coalescent.TreeIntervals;
import beast.util.TreeParser;
import junit.framework.TestCase;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class PreferentialBayesianSkylineTest extends TestCase {

    @Test
    public void test1 () {

        Tree tree = new TreeParser("((D4Philip56:2.0,(D4Philip64:3.0,D4Philip84:23.0):7.0):10.0,(D4SLanka78:19.0,(D4Thai78:5.0,D4Thai84:11.0):14.0):15.0);",false);
        TreeIntervals intervals = new TreeIntervals(tree);

        PreferentialBayesianSkyline skyline = new PreferentialBayesianSkyline();
        skyline.initByName("popSizes", "1.0 2.0 3.0", "groupSizes", "4 4 3", "treeIntervals", intervals, "samplingIntensity", "0.4");


        //double logP = skyline.calculateLogP();

        //System.out.println(logP);

    }


    /**********************************************************************/
    /* Comparisons to previous implmentation                              */
    /**********************************************************************/


    @Test
    public void testPSBSPTreeLikelihood1 () {

        System.out.println("Preferential sampling Bayesian Skyline compared to old Preferential Bayesian Skyline implementation: " +
                           "Tree with homochronous sampling.");

        Tree tree = new TreeParser("((D4Philip56:30.0,(D4Philip64:23.0,D4Philip84:23.0):7.0):10.0,(D4SLanka78:25.0,(D4Thai78:11.0,D4Thai84:11.0):14.0):15.0);",false);
        TreeIntervals intervals = new TreeIntervals(tree);

        GeneralBayesianSkyline skyline1 = new PreferentialSamplingBayesianSkyline();
        skyline1.initByName("popSizes", "1.0 2.0", "samplingIntensity", "2.0", "groupSizes", "8 3", "treeIntervals", intervals);

        PreferentialBayesianSkyline skyline2 = new PreferentialBayesianSkyline();
        skyline2.initByName("popSizes", "1.0 2.0", "samplingIntensity", "2.0", "groupSizes", "8 3", "treeIntervals", intervals);

        double logP1 = skyline1.calculateLogP();
        double logP2 = skyline2.calculateLogP();

        assertEquals(logP2, logP1);
        //System.out.println(logP1+"\t"+logP2);
    }


    @Test
    public void testPSBSPTreeLikelihood2 () {

        System.out.println("Preferential sampling Bayesian Skyline compared to old Preferential Bayesian Skyline implementation: " +
                "Tree with all sampling and coalescent times at unique times");

        Tree tree = new TreeParser("((((D4Mexico84:5.0,D4ElSal94:15.0):1.0,D4PRico86:8.0):1.0,D4Tahiti79:2.0):5.0,D4Indon77:5.0);",false);
        TreeIntervals intervals = new TreeIntervals(tree);

        GeneralBayesianSkyline skyline1 = new PreferentialSamplingBayesianSkyline();
        skyline1.initByName("popSizes", "1.0 2.0", "samplingIntensity", "2.0", "groupSizes", "4 5", "treeIntervals", intervals);

        PreferentialBayesianSkyline skyline2 = new PreferentialBayesianSkyline();
        skyline2.initByName("popSizes", "1.0 2.0", "samplingIntensity", "2.0", "groupSizes", "4 5", "treeIntervals", intervals);

        double logP1 = skyline1.calculateLogP();
        double logP2 = skyline2.calculateLogP();

        assertEquals(logP2, logP1);
        //System.out.println(logP1+"\t"+logP2);

    }


    @Test
    public void testBSPTreeLikelihood3 () {

        System.out.println("Preferential sampling Bayesian Skyline compared to old Preferential Bayesian Skyline implementation: " +
                "Tree with all coalescent times at unique times, sampling times not unique");

        Tree tree = new TreeParser("((D4Philip56:2.0,(D4Philip64:3.0,D4Philip84:23.0):7.0):10.0,(D4SLanka78:19.0,(D4Thai78:5.0,D4Thai84:11.0):14.0):15.0);",false);
        TreeIntervals intervals = new TreeIntervals(tree);

        GeneralBayesianSkyline skyline1 = new PreferentialSamplingBayesianSkyline();
        skyline1.initByName("popSizes", "1.0 2.0", "samplingIntensity", "2.0", "groupSizes", "8 3", "treeIntervals", intervals);

        PreferentialBayesianSkyline skyline2 = new PreferentialBayesianSkyline();
        skyline2.initByName("popSizes", "1.0 2.0", "samplingIntensity", "2.0", "groupSizes", "8 3", "treeIntervals", intervals);

        double logP1 = skyline1.calculateLogP();
        double logP2 = skyline2.calculateLogP();

        assertEquals(logP2, logP1);
        //System.out.println(logP1+"\t"+logP2);

    }


    @Test
    public void testBSPTreeLikelihood4 () {

        System.out.println("Preferential sampling Bayesian Skyline compared to old Preferential Bayesian Skyline implementation: " +
                "Tree with both coalescent and sampling times at non-unique times.");

        Tree tree = new TreeParser("((D4Philip56:2.0,(D4Philip64:3.0,D4Philip84:23.0):7.0):10.0,(D4SLanka78:17.0,(D4Thai78:5.0,D4Thai84:11.0):12.0):17.0);",false);
        TreeIntervals intervals = new TreeIntervals(tree);

        GeneralBayesianSkyline skyline1 = new PreferentialSamplingBayesianSkyline();
        skyline1.initByName("popSizes", "1.0 2.0", "samplingIntensity", "2.0", "groupSizes", "8 3", "treeIntervals", intervals);

        PreferentialBayesianSkyline skyline2 = new PreferentialBayesianSkyline();
        skyline2.initByName("popSizes", "1.0 2.0", "samplingIntensity", "2.0", "groupSizes", "8 3", "treeIntervals", intervals);

        double logP1 = skyline1.calculateLogP();
        double logP2 = skyline2.calculateLogP();

        assertEquals(logP2, logP1);
        //System.out.println(logP1+"\t"+logP2);

    }



}
