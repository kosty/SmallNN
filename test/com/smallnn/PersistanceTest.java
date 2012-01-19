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

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

import javax.vecmath.GMatrix;

import org.junit.Assert;
import org.junit.Test;

import com.smallnn.input.TrainDataUtil;
import com.smallnn.input.TrainDataUtil.Data;

public class PersistanceTest {

    @Test
    public void testPersistance() throws Exception {
        Data d = TrainDataUtil.blackWhiteData();
        NeuralNetwork sln0 = new SingleLayerNetwork(30, d.y.getNumCol(), d.x.getNumCol());
        
        sln0.train(d.x, d.y, 1, 0.6);
        
        GMatrix[] persisted = sln0.getConfig();

        
        ObjectOutputStream oos = null;
        try {
            oos = new ObjectOutputStream(new FileOutputStream("./configs/nhln.dat"));
            oos.writeObject(persisted);
        } finally {
            if (oos != null)
                oos.close();
        }
        
        GMatrix[] restored = null;
        ObjectInputStream ois = null;
        try{
            ois = new ObjectInputStream(new FileInputStream("./configs/nhln.dat"));
            restored = (GMatrix[]) ois.readObject();
        } finally {
            if (ois != null)
                ois.close();
        }
        for (int l = 0; l<persisted.length;l++){
            Assert.assertTrue(persisted[l].equals(restored[l]));
        }
        
        SingleLayerNetwork sln1 = new SingleLayerNetwork(restored[0], restored[1], restored[2], restored[3]);
        
    }

}
