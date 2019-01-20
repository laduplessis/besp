package bsp.distributions;

import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;

public class SimpleBSPDebugger extends SimpleBSP {

    final public Input<IntegerParameter> storesInput =
            new Input<>("stores","Number of calls to store()");

    final public Input<IntegerParameter> restoresInput =
            new Input<>("restores","Number of calls to restore()");

    final public Input<IntegerParameter> requiresRecalculationInput =
            new Input<>("requiresRecalculation","Number of calls to requiresRecalculation()");

    final public Input<IntegerParameter> updatesInput =
            new Input<>("updates","Number of calls to updateArrays()");


    int stores, restores, updates, recalcs;

    @Override
    public void initAndValidate() {
        stores = 0;
        restores = 0;
        updates = 0;
        recalcs = 0;

        storesInput.get().setValue(0,0);
        restoresInput.get().setValue(0,0);
        updatesInput.get().setValue(0,0);
        requiresRecalculationInput.get().setValue(0,0);

        super.initAndValidate();
    }



    @Override
    protected void updateArrays() {
        //System.out.println("Update!");
        updatesInput.get().setValue(0, ++updates);
        super.updateArrays();
    }

    @Override
    protected boolean requiresRecalculation() {
        requiresRecalculationInput.get().setValue(0, ++recalcs);
        return (super.requiresRecalculation());
    }

    @Override
    public void store() {
        storesInput.get().setValue(0, ++stores);
        super.store();
    }

    @Override
    public void restore() {
        restoresInput.get().setValue(0, ++restores);
        super.restore();
    }

}
