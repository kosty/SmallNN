package com.smallnn;

import javax.vecmath.GMatrix;

public interface NeuralNetwork {

    Double[] train(GMatrix x, GMatrix y, double lambda, double learningRate) throws Exception;

    void activate(GMatrix x);

    GMatrix getOutput();
    
    GMatrix[] getConfig();

}