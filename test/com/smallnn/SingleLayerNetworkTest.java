package com.smallnn;

import javax.vecmath.GMatrix;

import org.junit.Assert;
import org.junit.Test;

public class SingleLayerNetworkTest {

    @Test
    public void testActivate(){
        GMatrix theta1 = new GMatrix(1, 2, new double[]{0., 1.});
        GMatrix theta2 = new GMatrix(1, 2, new double[]{0., 1.});
        GMatrix mu     = new GMatrix(1, 2, new double[]{0., 0.});
        GMatrix sigma  = new GMatrix(1, 2, new double[]{1., 1.});
        SingleLayerNetwork sln = new SingleLayerNetwork(mu, sigma, theta1, theta2);
        
        GMatrix x = new GMatrix(4, 1, new double[]{1., 0., -1., 0.5});
//        GMatrix y = new GMatrix(4, 1, new double[]{1., 0.,  0., 1.});
        
        sln.activate(x);
        
        Assert.assertTrue(sln.A2.equals(new GMatrix(4,2, new double[]{1.0, 0.7310585786300049, 1.0, 0.5,  1.0, 0.2689414213699951, 1.0, 0.6224593312018546})));
        Assert.assertTrue(sln.Z3.equals(new GMatrix(4,1, new double[]{0.7310585786300049, 0.5, 0.2689414213699951, 0.6224593312018546})));
        Assert.assertTrue(sln.A3.equals(new GMatrix(4,1, new double[]{0.6750375273768237, 0.6224593312018546, 0.5668330070205946, 0.6507776782147005})));
        
        System.out.println("Success");
    }

}
