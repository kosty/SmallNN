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

import javax.vecmath.GMatrix;

public interface NeuralNetwork {

    /**
     * 
     * @param x
     * @param y
     * @param lambda
     * @param learningRate
     * @return
     * @throws Exception
     */
    Double[] train(GMatrix x, GMatrix y, double lambda, double learningRate) throws Exception;
    
    
    Double[] train(GMatrix x, GMatrix y, double precission, int maxIterations) throws Exception;

    /**
     * 
     * @param x
     * @return
     */
    GMatrix activate(GMatrix x);

    /**
     * 
     * @return
     */
    GMatrix[] getConfig();

}