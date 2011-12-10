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
