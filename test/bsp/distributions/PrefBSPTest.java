package bsp.distributions;

import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.TreeIntervals;
import beast.util.TreeParser;
import bsp.deprecated.EpochSamplingBayesianSkyline;
import junit.framework.TestCase;
import org.junit.Test;

public class PrefBSPTest extends TestCase {

    @Test
    public void test1 () {

        // Tree with 6 samples, 2 sampling events at time 0 and at time 6, all other samples and coalescences at unique times
        // Tree tree = new TreeParser("((D4Philip56:2.0,(D4Philip64:3.0,D4Philip84:23.0):7.0):10.0,(D4SLanka78:19.0,(D4Thai78:5.0,D4Thai84:11.0):14.0):15.0);",false);

        // Tree with 6 samples, 2 sampling events at time 0 and at time 6 coalescences at time 23
        // Tree tree = new TreeParser("((D4Philip56:2.0,(D4Philip64:3.0,D4Philip84:23.0):7.0):10.0,(D4SLanka78:17.0,(D4Thai78:5.0,D4Thai84:11.0):12.0):17.0);",false);

        // Tree with 7 samples, polytomy at time 11
        Tree tree = new TreeParser("((D4Philip56:2.0,(D4Philip64:3.0,D4Philip84:23.0):7.0):10.0,(D4SLanka78:19.0,(D4Thai78:5.0,D4Thai84:11.0):14.0):15.0);",false);
        TreeIntervals intervals = new TreeIntervals(tree);
        double [] times = new double [] {9,10,15,20};



        BSP skyline = new PrefBSP();
        skyline.initByName("treeIntervals", intervals, "popSizes", "1.0 2.0 3.0", "popSizeGroupSizes", "2 4 5",
                           "samplingIntensity", "2 3", "samplingIntensityGroupSizes", "2 4");

        System.out.println("Intervals: "+intervals.getIntervalCount());
        System.out.println("Coalescences: "+intervals.getSampleCount());

        //skyline.collectSegmentTimes();

        System.out.println(skyline);
        System.out.println(skyline.calculateLogP());

    }



    @Test
    public void testPSBSPTreeLikelihood1 () {

        System.out.println("Preferential sampling Bayesian Skyline compared to old Preferential Bayesian Skyline implementation: " +
                "Tree with homochronous sampling.");

        Tree tree = new TreeParser("((D4Philip56:30.0,(D4Philip64:23.0,D4Philip84:23.0):7.0):10.0,(D4SLanka78:25.0,(D4Thai78:11.0,D4Thai84:11.0):14.0):15.0);",false);
        TreeIntervals intervals = new TreeIntervals(tree);

        GeneralBayesianSkyline skyline1 = new PreferentialSamplingBayesianSkyline();
        skyline1.initByName("popSizes", "1.0 2.0", "samplingIntensity", "2.0", "groupSizes", "8 3", "treeIntervals", intervals);

        BSP skyline2 = new PrefBSP();
        skyline2.initByName("popSizes", "1.0 2.0", "samplingIntensity", "2.0", "popSizeGroupSizes", "8 3", "treeIntervals", intervals);

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

        BSP skyline2 = new PrefBSP();
        skyline2.initByName("popSizes", "1.0 2.0", "samplingIntensity", "2.0", "popSizeGroupSizes", "4 5", "treeIntervals", intervals);


        //PreferentialBayesianSkyline skyline2 = new PreferentialBayesianSkyline();
        //skyline2.initByName("popSizes", "1.0 2.0", "samplingIntensity", "2.0", "groupSizes", "4 5", "treeIntervals", intervals);

        System.out.println(skyline2);

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

        BSP skyline2 = new PrefBSP();
        skyline2.initByName("popSizes", "1.0 2.0", "samplingIntensity", "2.0", "popSizeGroupSizes", "8 3", "treeIntervals", intervals);

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

        BSP skyline2 = new PrefBSP();
        skyline2.initByName("popSizes", "1.0 2.0", "samplingIntensity", "2.0", "popSizeGroupSizes", "8 3", "treeIntervals", intervals);

        double logP1 = skyline1.calculateLogP();
        double logP2 = skyline2.calculateLogP();

        assertEquals(logP2, logP1);
        //System.out.println(logP1+"\t"+logP2);

    }


    @Test
    public void testPSBSPTreeEpochTest () {

        System.out.println("Preferential sampling Bayesian Skyline compared to old Preferential Bayesian Skyline implementation: " +
                "Tree with all sampling and coalescent times at unique times");

        Tree tree = new TreeParser("((((D4Mexico84:5.0,D4ElSal94:15.0):1.0,D4PRico86:8.0):1.0,D4Tahiti79:2.0):5.0,D4Indon77:5.0);",false);
        TreeIntervals intervals = new TreeIntervals(tree);

        BSP skyline2 = new PrefBSP();
        skyline2.initByName("popSizes", "1.0 2.0", "samplingIntensity", "1.0 2.0 3.0", "popSizeGroupSizes", "4 5",
                            "samplingIntensityGroupSizes", "1 1 1", "treeIntervals", intervals);


        //PreferentialBayesianSkyline skyline2 = new PreferentialBayesianSkyline();
        //skyline2.initByName("popSizes", "1.0 2.0", "samplingIntensity", "2.0", "groupSizes", "4 5", "treeIntervals", intervals);

        System.out.println(skyline2);

        double logP2 = skyline2.calculateLogP();

        System.out.println(logP2);

    }

}