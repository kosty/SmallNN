package com.smallnn;

import javax.vecmath.GMatrix;

public interface NeuralNetwork {

    Double[] train(GMatrix x, GMatrix y, double lambda, double learningRate) throws Exception;

    GMatrix activate(GMatrix x);

    GMatrix[] getConfig();

}