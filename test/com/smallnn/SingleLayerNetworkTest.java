/**
 *  This file is part of SmallNN, a small neural network implementation
 *  Copyright (C) 2011, 2012 Arsen Kostenko <arsen.kostenko@gmail.com>
 *     
 *  SmallNN is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  SmallNN is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with SmallNN.  If not, see <http://www.gnu.org/licenses/>.
 */
package com.smallnn;

import static org.junit.Assert.assertTrue;

import javax.vecmath.GMatrix;

import org.junit.Test;

public class SingleLayerNetworkTest {

    @Test
    public void testActivate(){
        GMatrix theta1 = new GMatrix(1, 2, new double[]{0., 1.});
        GMatrix theta2 = new GMatrix(1, 2, new double[]{0., 1.});
        GMatrix mu     = new GMatrix(1, 1, new double[]{0.});
        GMatrix sigma  = new GMatrix(1, 1, new double[]{1.});
        
        System.out.println(mu.getNumCol()+" vs "+(theta1.getNumCol()-1));
        
        SingleLayerNetwork sln = new SingleLayerNetwork(mu, sigma, theta1, theta2);
        
        GMatrix x = new GMatrix(4, 1, new double[]{1., 0., -1., 0.5});
//        GMatrix y = new GMatrix(4, 1, new double[]{1., 0.,  0., 1.});
        
        sln.activate(x);
        
        assertTrue(sln.A2.equals(new GMatrix(4,2, new double[]{1.0, 0.7310585786300049, 1.0, 0.5,  1.0, 0.2689414213699951, 1.0, 0.6224593312018546})));
        assertTrue(sln.Z3.equals(new GMatrix(4,1, new double[]{0.7310585786300049, 0.5, 0.2689414213699951, 0.6224593312018546})));
        assertTrue(sln.A3.equals(new GMatrix(4,1, new double[]{0.6750375273768237, 0.6224593312018546, 0.5668330070205946, 0.6507776782147005})));
        
        System.out.println("Success");
    }

}
