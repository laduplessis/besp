package bsp.distributions;

import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.TreeIntervals;
import beast.util.TreeParser;
import junit.framework.TestCase;
import org.junit.Test;
import test.beast.BEASTTestCase;

/**
 * Unit tests for bsp.distributions.BESP.java
 *
 *     1. Homochronous tree - Done
 *
 * - Compare likelihood to results calculated by hand for:
 *      2. Heterochronous tree (all times unique)
 *      3. Heterochronous tree (non-unique sampling times)
 *      4. Heterochronous tree (non-unique sampling and coalescent times)
 *      5. Heterochronous tree (sampling event between coalescent event between segments)
 *
 * Polytomy?
 */

public class BESPTest extends TestCase {


    /*********************/
    /* Comparison to BSP */
    /*********************/

    @Test
    public void testBESPTreeLikelihood1 () {

        System.out.println("BESP compared to BSP: " +
                "Tree with homochronous sampling (likelihood should be identical).");

        Tree tree = new TreeParser("((D4Philip56:30.0,(D4Philip64:23.0,D4Philip84:23.0):7.0):10.0,(D4SLanka78:25.0,(D4Thai78:11.0,D4Thai84:11.0):14.0):15.0);",false);
        TreeIntervals intervals = new TreeIntervals(tree);

        BSP skyline1 = new BSP();
        skyline1.initByName("popSizes", "1.0 2.0", "popSizeGroupSizes", "2 3", "treeIntervals", intervals);

        // Need to add all sampling events (6) to the first popSize segment
        BSP skyline2 = new BESP();
        skyline2.initByName("popSizes", "1.0 2.0", "popSizeGroupSizes", "8 3", "samplingIntensity", "1.0", "treeIntervals", intervals);


        double logP1 = skyline1.calculateLogP();
        double logP2 = skyline2.calculateLogP();

        assertEquals(logP2, logP1);
        System.out.println(logP1+"\t"+logP2);
    }


    /************************************************/
    /* Comparisons to likelihood calculated by hand */
    /************************************************/


    @Test
    public void testBESPTreeLikelihood2 () {

        System.out.println("BSP compared to old Bayesian Skyline implementation: "+
                "Tree with all sampling and coalescent times at unique times");

        Tree tree = new TreeParser("((((D4Mexico84:5.0,D4ElSal94:15.0):1.0,D4PRico86:8.0):1.0,D4Tahiti79:2.0):5.0,D4Indon77:5.0);",false);
        TreeIntervals intervals = new TreeIntervals(tree);

        BSP skyline1 = new BESP();
        skyline1.initByName("popSizes", "1.0 2.0", "samplingIntensity", "2.0", "popSizeGroupSizes", "6 3", "treeIntervals", intervals);

        //System.out.println(skyline1);

        double logP1 = skyline1.calculateLogP();

        assertEquals(-56.2274112777602, logP1, BEASTTestCase.PRECISION);
        // System.out.println(logP1);
    }

    @Test
    public void testBESPTreeLikelihood3 () {

        System.out.println("BSP compared to old Bayesian Skyline implementation: "+
                "Tree with all coalescent times at unique times, sampling times not unique");

        Tree tree = new TreeParser("((D4Philip56:2.0,(D4Philip64:3.0,D4Philip84:23.0):7.0):10.0,(D4SLanka78:19.0,(D4Thai78:5.0,D4Thai84:11.0):14.0):15.0);",false);
        TreeIntervals intervals = new TreeIntervals(tree);

        BSP skyline1 = new BESP();
        skyline1.initByName("popSizes", "2.0 3.0", "samplingIntensity", "1.0 2.0 3.0",
                            "popSizeGroupSizes", "6 5", "samplingIntensityGroupSizes", "3 2 1", "treeIntervals", intervals);

        System.out.println(skyline1);

        double logP1 = skyline1.calculateLogP();

        assertEquals(-183.8716748273099, logP1, BEASTTestCase.PRECISION);
        System.out.println(logP1);
    }

    @Test
    public void testBESPTreeLikelihood4 () {

        System.out.println("BSP compared to old Bayesian Skyline implementation: "+
                "Tree with both coalescent and sampling times at non-unique times.");

        Tree tree = new TreeParser("((D4Philip56:2.0,(D4Philip64:3.0,D4Philip84:23.0):7.0):10.0,(D4SLanka78:17.0,(D4Thai78:5.0,D4Thai84:11.0):12.0):17.0);",false);
        TreeIntervals intervals = new TreeIntervals(tree);

        BSP skyline1 = new BESP();
        skyline1.initByName("popSizes", "3.0 2.0 1.0", "samplingIntensity", "2.0 3.0",
                "popSizeGroupSizes", "3 4 4", "samplingIntensityGroupSizes", "4 2", "treeIntervals", intervals);

        System.out.println(skyline1);

        double logP1 = skyline1.calculateLogP();

        assertEquals(-205.234349834419, logP1, BEASTTestCase.PRECISION);
    }

    @Test
    public void testBESPTreeLikelihood5 () {

        System.out.println("BSP compared to old Bayesian Skyline implementation: "+
                "Tree with a sampling event between the coalescent events that define a group boundary");

        Tree tree = new TreeParser("((D4Philip56:2.0,(D4Philip64:3.0,D4Philip84:23.0):7.0):10.0,(D4SLanka78:26.0,(D4Thai78:5.0,D4Thai84:11.0):21.0):8.0);",false);
        TreeIntervals intervals = new TreeIntervals(tree);

        BSP skyline1 = new BESP();
        skyline1.initByName("popSizes", "2.0 3.0 4.0", "samplingIntensity", "2.0 3.0 4.0",
                "popSizeGroupSizes", "5 2 4", "samplingIntensityGroupSizes", "3 2 1", "treeIntervals", intervals);

        System.out.println(skyline1);

        double logP1 = skyline1.calculateLogP();

        assertEquals(-289.280186700424, logP1, BEASTTestCase.PRECISION);
    }



}