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

import java.util.Arrays;

public class PushOutQueue {
    
    double[] queue;
    
    public PushOutQueue(int size, double v){
        this.queue = new double[size];
        Arrays.fill(this.queue, v);
    }
    
    public void add(double d){
        pushOut();
        this.queue[this.queue.length-1] = d;
    }

    private void pushOut() {
        for(int i=0; i < this.queue.length-1; i++)
            this.queue[i] = this.queue[i+1];
    }
    
    public double peek(int idx){
        return this.queue[idx];
    }

}
