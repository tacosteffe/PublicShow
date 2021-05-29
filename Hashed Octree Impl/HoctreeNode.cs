using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public enum NODE_TYPE
{
    NONE,
    INTERNAL,
    PSUEDO,
    LEAF
}

public class HoctreeNode
{
    public NODE_TYPE Type;
    public uint LocCode; //LocationCode, contains min and depth/size of node

    public VoxelData Data;


    public HoctreeNode()
    {
        LocCode = 0;
        Data.Qef = new QefData();
    }

    public uint GetID
    {
        get { return LocCode & 7; }
    }

}

//Data partaining the information that is needed for dual contouring and voxel modifications
public struct VoxelData
{
    //Direct voxel attributes.
    //Might be able to generate these on a need to know basis

    //The voxel type partaining the information about the voxel and draw information
    public COLOR VoxelType;
    //The voxel density, needed for deformation purposes
    public float VoxelDensity;



    //DC needed objects
    public Vector3 Position; //Calculated position of the vertex?
    public Vector3 Normal; //The normal of the vertex.. Might be able to remove this and the AvgNormal
    public Vector3 AverageNormal;
    public int Index; //Index used for generating mesh, might be able to remove this
    public int Corners; //Corners of the Voxel Cube, same here might be able to remove this
    public QefData Qef; //QEF data, used for DC contouring and simplification.. lots of data i wish to get rid off
}

public enum COLOR
{
    NONE = 0,
    GREEN,
    BROWN,
    GREY
}

public static class VoxelExtensions
{
   //Testing for vertex colors.
    public static Color GetColor(this COLOR vox)
    {

        switch(vox)
        {
            case COLOR.NONE:
                return Color.magenta;

            case COLOR.GREEN:
                return Color.green;

            case COLOR.BROWN:
                return new Color(165f / 255f, 42f / 255f, 42f / 255f); //BROWN

            case COLOR.GREY:
                return Color.grey;

            default:
                return Color.magenta;
        }
    }

}