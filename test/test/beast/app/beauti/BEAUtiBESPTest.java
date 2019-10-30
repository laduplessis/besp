package test.beast.app.beauti;

import org.fest.swing.fixture.JTabbedPaneFixture;
import org.junit.Test;

import java.io.File;

/* Run from beast2 module path (examples/nexus/Flu.nex needs to be in path) */
public class BEAUtiBESPTest extends BeautiBase {

    @Test
    public void beautiBSPTest() throws Exception {
        warning("Load Flu.nex");
        importAlignment("examples/nexus", new File("Flu.nex"));

        JTabbedPaneFixture f = beautiFrame.tabbedPane();
        f.selectTab("Priors");

        setAndCheckBSPTreePrior(f);
        setAndCheckExpGrowthTreePrior(f);
        setAndCheckBSPTreePrior(f);

        makeSureXMLParses();
    }

    @Test
    public void beautiBESPTest() throws Exception {
        warning("Load Flu.nex");
        importAlignment("examples/nexus", new File("Flu.nex"));

        JTabbedPaneFixture f = beautiFrame.tabbedPane();
        f.selectTab("Priors");

        setAndCheckBESPTreePrior(f);
        setAndCheckExpGrowthTreePrior(f);
        setAndCheckBESPTreePrior(f);
        setAndCheckBSPTreePrior(f);
        setAndCheckBESPTreePrior(f);

        makeSureXMLParses();
    }



    /* Helper functions to set and check tree priors */

    private void setAndCheckBSPTreePrior(JTabbedPaneFixture f) throws Exception {

        warning("Change to Coalescent - BSP");
        beautiFrame.comboBox("TreeDistribution").selectItem("Coalescent BSP");
        printBeautiState(f);

        assertStateEquals("Tree.t:Flu", "BSPPopSizes.t:Flu", "BSPGroupSizes.t:Flu");
        assertOperatorsEqual( "BSPTreeScaler.t:Flu",
                "BSPTreeRootScaler.t:Flu",
                "BSPUniformOperator.t:Flu",
                "BSPSubtreeSlide.t:Flu",
                "BSPNarrow.t:Flu",
                "BSPWide.t:Flu",
                "BSPWilsonBalding.t:Flu",
                "BSPPopSizesScaler.t:Flu",
                "BSPGroupSizesDelta.t:Flu");
        assertPriorsEqual("BSP.t:Flu", "BSPMarkovChainedPopSizes.t:Flu");
        assertTraceLogEqual(  "posterior",
                "likelihood",
                "prior",
                "treeLikelihood.Flu",
                "TreeHeight.t:Flu",
                "BSP.t:Flu",
                "BSPPopSizes.t:Flu",
                "BSPGroupSizes.t:Flu",
                "BSPPopSizeChangeTimes.t:Flu");

        makeSureXMLParses();

    }


    private void setAndCheckBESPTreePrior(JTabbedPaneFixture f) throws Exception {

        warning("Change to Coalescent - BESP");
        beautiFrame.comboBox("TreeDistribution").selectItem("Coalescent BESP");
        printBeautiState(f);

        assertStateEquals("Tree.t:Flu",
                "BESPPopSizes.t:Flu",
                "BESPGroupSizes.t:Flu",
                "BESPSamplingIntensity.t:Flu");
        assertOperatorsEqual( "BESPTreeScaler.t:Flu",
                "BESPTreeRootScaler.t:Flu",
                "BESPUniformOperator.t:Flu",
                "BESPSubtreeSlide.t:Flu",
                "BESPNarrow.t:Flu",
                "BESPWide.t:Flu",
                "BESPWilsonBalding.t:Flu",
                "BESPPopSizesScaler.t:Flu",
                "BESPGroupSizesDelta.t:Flu",
                "BESPSamplingIntensityScaler.t:Flu");
        assertPriorsEqual("BESP.t:Flu",
                "BESPMarkovChainedPopSizes.t:Flu",
                "BESPSamplingIntensityPrior.t:Flu");
        assertTraceLogEqual(  "posterior",
                "likelihood",
                "prior",
                "treeLikelihood.Flu",
                "TreeHeight.t:Flu",
                "BESP.t:Flu",
                "BESPPopSizes.t:Flu",
                "BESPGroupSizes.t:Flu",
                "BESPSamplingIntensity.t:Flu",
                "BESPPopSizeChangeTimes.t:Flu",
                "BESPSamplingIntensityChangeTimes.t:Flu");

        makeSureXMLParses();

    }


    private void setAndCheckExpGrowthTreePrior(JTabbedPaneFixture f) throws Exception {

        warning("Change to Coalescent - exponential population");
        beautiFrame.comboBox("TreeDistribution").selectItem("Coalescent Exponential Population");
        printBeautiState(f);

        assertStateEquals("Tree.t:Flu", "ePopSize.t:Flu", "growthRate.t:Flu");
        assertOperatorsEqual( "CoalescentExponentialTreeScaler.t:Flu",
                "CoalescentExponentialTreeRootScaler.t:Flu",
                "CoalescentExponentialUniformOperator.t:Flu",
                "CoalescentExponentialSubtreeSlide.t:Flu",
                "CoalescentExponentialNarrow.t:Flu",
                "CoalescentExponentialWide.t:Flu",
                "CoalescentExponentialWilsonBalding.t:Flu",
                "ePopSizeScaler.t:Flu",
                "GrowthRateRandomWalk.t:Flu");
        assertPriorsEqual("CoalescentExponential.t:Flu", "ePopSizePrior.t:Flu", "GrowthRatePrior.t:Flu");
        assertTraceLogEqual(  "posterior",
                "likelihood",
                "prior",
                "treeLikelihood.Flu",
                "TreeHeight.t:Flu",
                "CoalescentExponential.t:Flu",
                "ePopSize.t:Flu",
                "growthRate.t:Flu");

        makeSureXMLParses();
    }




}
