using System.Collections;
using System.Collections.Generic;
using UnityEngine;


//NOTES:
/*
    Most of the tables and Octree handling are taken from NGildea's blog.
    I claim no credit to those parts.
    
    Source:
    https://ngildea.blogspot.com/2014/11/implementing-dual-contouring.html
    https://ngildea.blogspot.com/2014/09/dual-contouring-chunked-terrain.html

    C# Port: -Tuckbone
    https://github.com/tuckbone/DualContouringCSharp

    Shabby Hashed octree impl:
    Yours truly.. <- remove
 */

//Hashed octree impl.
//Think of this octree as a chunk.
//REDUDANT: Skip octree implementation at whole.
/*
    Annoying handling of voxels, in performance sure it's good enough but editing is less nice when working at corners of chunks.
    
    Pros/cons for future me to worry about

    PROS:
    1: Good partitioning when done at full
    2: Easier handling of merge of several chunks, but at the same time might be worse.. dunno have to test

    CONS:
    1: turn Nodes to nullable structs, overhead from class takes a bit when working in 1000s
    2: Needs simplification when getting a voxel value from coords.
    3: Memory trade of is time taken to get the voxel in question. +/- depends on software and use i suppose
 */
public class Hoctree
{

    #region Octree Implementation

    //Pretty much the lifted version of NGildea's octree implementation ported to C# by 



    private Vector3 Position; //Position will be the min (0,0,0) of the octree
    //Hoctree with UINT can handle a maximum of 11 in depth
    private short Size = 16; //Chunk size
    private Hashtable Nodes; //Child nodes

    private GameObject MeshObject = null;


    public Hoctree(Vector3 pos, short size)
    {
        Position = pos;
        Size = size;
        Nodes = new Hashtable();
    }

    public HoctreeNode GetNode(uint locCode)
    {
        if(Nodes.ContainsKey(locCode))
            return Nodes[locCode] as HoctreeNode;
        return null;
    }


    public void ConstructOctree()
    {
        HoctreeNode parent = new HoctreeNode();
        parent.LocCode = 1;
        parent.Type = NODE_TYPE.INTERNAL;
        Nodes.Add(parent.LocCode, parent);

        for(uint i = 0; i < 8; ++i)
        {
            //First time nodes work a bit different in LocCode
            //They will get assigned a 1 in the in the beggining of their LocCode
            ConstructNode(i, parent.LocCode, Size / 2);
        }

        //TODO not working
        //Nodes[parent.LocCode] = SimplifyOctree(parent, .1f);

        List<MeshVertex> vertexBuffer = new List<MeshVertex>();
        List<int> indexBuffer = new List<int>();
        GenerateMeshFromOctree(parent, vertexBuffer, indexBuffer);
    }

    public void DeleteHoctree()
    {
        GameObject.Destroy(MeshObject);
    }

    public HoctreeNode ConstructNode(uint id, uint parentId, int nSize)
    {
        HoctreeNode node = new HoctreeNode();
        node.LocCode = (parentId << 3); //Push parent ID back 3 bits
        node.LocCode += id; //Location code, follows CHILD_OFFSETS in ending bitwise

        if (nSize > 1)
        {
            node.Type = NODE_TYPE.INTERNAL;
            for (uint i = 0; i < 8; ++i)
            {
                ConstructNode(i, node.LocCode, nSize / 2);
            }
        }
        else if(nSize == 1)
        {
            ConstructLeaf(node);
        }

        if(node != null)
            Nodes.Add(node.LocCode, node);

        return node;
    }


    public HoctreeNode ConstructLeaf(HoctreeNode leaf)
    {
        Vector3 min = GetMinFromNode(leaf);//leaf.Min;
        int corners = 0;
        for (int i = 0; i < 8; i++)
        {
            Vector3 cornerPos = min + CHILD_MIN_OFFSETS[i];
            float density = glm.Density_Func(cornerPos);
            int material = density < 0.0f ? MATERIAL_SOLID : MATERIAL_AIR;
            corners |= (material << i);
        }

        if (corners == 0 || corners == 255)
        {
            // voxel is full inside or outside the volume
            //delete leaf
            //setting as null isn't required by the GC in C#... but its in the original, so why not!
            leaf = null; //Explicit deref
            return null;
        }

        // otherwise the voxel contains the surface, so find the edge intersections
        const int MAX_CROSSINGS = 6;
        int edgeCount = 0;
        Vector3 averageNormal = Vector3.zero;
        QefSolver qef = new QefSolver();

        for (int i = 0; i < 12 && edgeCount < MAX_CROSSINGS; i++)
        {
            int c1 = edgevmap[i][0];
            int c2 = edgevmap[i][1];

            int m1 = (corners >> c1) & 1;
            int m2 = (corners >> c2) & 1;

            if ((m1 == MATERIAL_AIR && m2 == MATERIAL_AIR) || (m1 == MATERIAL_SOLID && m2 == MATERIAL_SOLID))
            {
                // no zero crossing on this edge
                continue;
            }

            Vector3 p1 = min + CHILD_MIN_OFFSETS[c1];
            Vector3 p2 = min + CHILD_MIN_OFFSETS[c2];
            Vector3 p = ApproximateZeroCrossingPosition(p1, p2);
            Vector3 n = CalculateSurfaceNormal(p);
            qef.add(p.x, p.y, p.z, n.x, n.y, n.z);

            averageNormal += n;

            edgeCount++;
        }

        Vector3 qefPosition = Vector3.zero;
        qef.solve(qefPosition, QEF_ERROR, QEF_SWEEPS, QEF_ERROR);

        leaf.Data.Corners = 0;
        leaf.Data.Index = -1;
        leaf.Data.Position = new Vector3(qefPosition.x, qefPosition.y, qefPosition.z);
        leaf.Data.Qef = qef.getData();

        int nSize = GetDepthAndSize(leaf)[1]; //leaf.Size;
        Vector3 max = new Vector3(min.x + nSize, min.y + nSize, min.z + nSize);
        if (leaf.Data.Position.x < min.x || leaf.Data.Position.x > max.x ||
            leaf.Data.Position.y < min.y || leaf.Data.Position.y > max.y ||
            leaf.Data.Position.z < min.z || leaf.Data.Position.z > max.z)
        {
            leaf.Data.Position = qef.getMassPoint();
        }

        leaf.Data.AverageNormal = Vector3.Normalize(averageNormal / (float)edgeCount);
        leaf.Data.Corners = corners;
        leaf.Data.VoxelType = Random.Range(0, 3) == 0 ? COLOR.GREEN : COLOR.BROWN;

        leaf.Type = NODE_TYPE.LEAF;

        return leaf;
    }


    public void DestroyOctree(HoctreeNode node)
    {
        if (node == null)
        {
            return;
        }

        HoctreeNode[] childs = GetChildsFromNode(node);
        for (int i = 0; i < 8; i++)
        {
            DestroyOctree(childs[i]);
        }

        Nodes.Remove(node.LocCode);
    }

    public void GenerateMeshFromOctree(HoctreeNode node, List<MeshVertex> vertexBuffer, List<int> indexBuffer)
    {
        if (node == null)
        {
            return;
        }

        vertexBuffer = new List<MeshVertex>();
        indexBuffer = new List<int>();

        GenerateVertexIndices(node, vertexBuffer);
        ContourCellProc(node, indexBuffer);

        MeshObject = new GameObject("Mesh");
        Mesh mesh;
        MeshFilter filter;
        MeshRenderer meshRenderer;

        meshRenderer = MeshObject.AddComponent<MeshRenderer>();
        filter = MeshObject.AddComponent<MeshFilter>();

        //NOTE: Fling in a material you want to use.
        //For this solution it needs a vertex color shader attached to it as well if you want colors
        //meshRenderer.sharedMaterial = 

        Vector3[] vertArray = new Vector3[vertexBuffer.Count];
        //Vector2[] uvs = new Vector2[vertexBuffer.Count];
        Color[] vertColors = new Color[vertexBuffer.Count];
        for (int i = 0; i < vertexBuffer.Count; i++)
        {
            vertArray[i] = vertexBuffer[i].xyz;
            //uvs[i] = new Vector2(vertexBuffer[i].xyz.x, vertexBuffer[i].xyz.z); //No need
            vertColors[i] = vertexBuffer[i].color;
        }

        Vector3[] normsArray = new Vector3[vertexBuffer.Count];
        for (int i = 0; i < vertexBuffer.Count; i++)
        {
            normsArray[i] = vertexBuffer[i].normal;
        }

        mesh = filter.mesh;
        mesh.vertices = vertArray;
        //mesh.uv = uvs;
        mesh.SetColors(vertColors);
        mesh.triangles = indexBuffer.ToArray();
        mesh.normals = normsArray;
        mesh.RecalculateBounds();
    }


    public Vector3 ApproximateZeroCrossingPosition(Vector3 p0, Vector3 p1)
    {
        // approximate the zero crossing by finding the min value along the edge
        float minValue = 100000f;
        float t = 0f;
        float currentT = 0f;
        const int steps = 8;
        const float increment = 1f / (float)steps;
        while (currentT <= 1.0f)
        {
            Vector3 p = p0 + ((p1 - p0) * currentT);
            float density = Mathf.Abs(glm.Density_Func(p));
            if (density < minValue)
            {
                minValue = density;
                t = currentT;
            }

            currentT += increment;
        }

        return p0 + ((p1 - p0) * t);
    }

    public Vector3 CalculateSurfaceNormal(Vector3 p)
    {
        float H = 0.001f;
        float dx = glm.Density_Func(p + new Vector3(H, 0.0f, 0.0f)) - glm.Density_Func(p - new Vector3(H, 0.0f, 0.0f));
        float dy = glm.Density_Func(p + new Vector3(0.0f, H, 0.0f)) - glm.Density_Func(p - new Vector3(0.0f, H, 0.0f));
        float dz = glm.Density_Func(p + new Vector3(0.0f, 0.0f, H)) - glm.Density_Func(p - new Vector3(0.0f, 0.0f, H));

        return new Vector3(dx, dy, dz).normalized;
    }


    public void GenerateVertexIndices(HoctreeNode node, List<MeshVertex> vertexBuffer)
    {
        if (node == null)
        {
            return;
        }

        if (node.Type != NODE_TYPE.LEAF)
        {
            uint[] childs = GetChildLocsFromNode(node);
            for (int i = 0; i < 8; i++)
            {
                GenerateVertexIndices(GetNode(childs[i]), vertexBuffer);
            }
        }

        if (node.Type != NODE_TYPE.INTERNAL)
        {
            node.Data.Index = vertexBuffer.Count;
            vertexBuffer.Add(new MeshVertex(node.Data.Position, node.Data.AverageNormal, node.Data.VoxelType.GetColor()));
        }
    }

    //Unless you want to create your own DC table treat this like number magic.
    public void ContourProcessEdge(HoctreeNode[] node, int dir, List<int> indexBuffer)
    {
        int minSize = 1000000;		// arbitrary big number
        int minIndex = 0;
        int[] indices = new int[4] { -1, -1, -1, -1 };
        bool flip = false;
        bool[] signChange = new bool[4] { false, false, false, false };

        for (int i = 0; i < 4; i++)
        {
            int edge = processEdgeMask[dir][i];
            int c1 = edgevmap[edge][0];
            int c2 = edgevmap[edge][1];

            int m1 = (node[i].Data.Corners >> c1) & 1;
            int m2 = (node[i].Data.Corners >> c2) & 1;

            int nSize = GetDepthAndSize(node[i])[1]; //node[i].Size;
            if (nSize < minSize)
            {
                minSize = nSize;
                minIndex = i;
                flip = m1 != MATERIAL_AIR;
            }

            indices[i] = node[i].Data.Index;

            signChange[i] = m1 != m2;
        }

        if (signChange[minIndex])
        {
            if (!flip)
            {
                indexBuffer.Add(indices[0]);
                indexBuffer.Add(indices[1]);
                indexBuffer.Add(indices[3]);

                indexBuffer.Add(indices[0]);
                indexBuffer.Add(indices[3]);
                indexBuffer.Add(indices[2]);
            }
            else
            {
                indexBuffer.Add(indices[0]);
                indexBuffer.Add(indices[3]);
                indexBuffer.Add(indices[1]);

                indexBuffer.Add(indices[0]);
                indexBuffer.Add(indices[2]);
                indexBuffer.Add(indices[3]);
            }

        }
    }

    //Unless you want to create your own DC table treat this like number magic.
    public void ContourEdgeProc(HoctreeNode[] node, int dir, List<int> indexBuffer)
    {
        if (node[0] == null || node[1] == null || node[2] == null || node[3] == null)
        {
            return;
        }

        if (node[0].Type != NODE_TYPE.INTERNAL &&
            node[1].Type != NODE_TYPE.INTERNAL &&
            node[2].Type != NODE_TYPE.INTERNAL &&
            node[3].Type != NODE_TYPE.INTERNAL)
        {
            ContourProcessEdge(node, dir, indexBuffer);
        }
        else
        {
            for (int i = 0; i < 2; i++)
            {
                HoctreeNode[] edgeNodes = new HoctreeNode[4];
                int[] c = new int[4]
            {
                edgeProcEdgeMask[dir][i][0],
                edgeProcEdgeMask[dir][i][1],
                edgeProcEdgeMask[dir][i][2],
                edgeProcEdgeMask[dir][i][3],
            };

                for (int j = 0; j < 4; j++)
                {
                    uint[] childs = GetChildLocsFromNode(node[j]);
                    if (node[j].Type == NODE_TYPE.LEAF || node[j].Type == NODE_TYPE.PSUEDO)
                    {
                        edgeNodes[j] = node[j];
                    }
                    else
                    {
                        edgeNodes[j] = GetNode(childs[c[j]]);
                    }
                }

                ContourEdgeProc(edgeNodes, edgeProcEdgeMask[dir][i][4], indexBuffer);
            }
        }
    }

    //Unless you want to create your own DC table treat this like number magic.
    public void ContourFaceProc(HoctreeNode[] node, int dir, List<int> indexBuffer)
    {
        if (node[0] == null || node[1] == null)
        {
            return;
        }

        if (node[0].Type == NODE_TYPE.INTERNAL ||
            node[1].Type == NODE_TYPE.INTERNAL)
        {
            for (int i = 0; i < 4; i++)
            {
                HoctreeNode[] faceNodes = new HoctreeNode[2];
                int[] c = new int[2]
                {
                    faceProcFaceMask[dir][i][0],
                    faceProcFaceMask[dir][i][1],
                };

                for (int j = 0; j < 2; j++)
                {
                    uint[] childs = GetChildLocsFromNode(node[j]);
                    if (node[j].Type != NODE_TYPE.INTERNAL)
                    {
                        faceNodes[j] = node[j];
                    }
                    else
                    {
                        faceNodes[j] = GetNode(childs[c[j]]);
                    }
                }

                ContourFaceProc(faceNodes, faceProcFaceMask[dir][i][2], indexBuffer);
            }

            int[][] orders = new int[2][]
            {
                new int[4]{ 0, 0, 1, 1 },
                new int[4]{ 0, 1, 0, 1 },
            };

            for (int i = 0; i < 4; i++)
            {
                HoctreeNode[] edgeNodes = new HoctreeNode[4];
                int[] c = new int[4]
                {
                    faceProcEdgeMask[dir][i][1],
                    faceProcEdgeMask[dir][i][2],
                    faceProcEdgeMask[dir][i][3],
                    faceProcEdgeMask[dir][i][4],
                };

                int[] order = orders[faceProcEdgeMask[dir][i][0]];
                for (int j = 0; j < 4; j++)
                {

                    if (node[order[j]].Type == NODE_TYPE.LEAF ||
                        node[order[j]].Type == NODE_TYPE.PSUEDO)
                    {
                        edgeNodes[j] = node[order[j]];
                    }
                    else
                    {
                        uint[] childs = GetChildLocsFromNode(node[order[j]]);
                        edgeNodes[j] = GetNode(childs[c[j]]);
                    }
                }

                ContourEdgeProc(edgeNodes, faceProcEdgeMask[dir][i][5], indexBuffer);
            }
        }
    }

    //Unless you want to create your own DC table treat this like number magic.
    public void ContourCellProc(HoctreeNode node, List<int> indexBuffer)
    {
        if (node == null)
        {
            return;
        }

        if (node.Type == NODE_TYPE.INTERNAL)
        {
            uint[] childs = GetChildLocsFromNode(node);

            for (int i = 0; i < 8; i++)
            {
                ContourCellProc(GetNode(childs[i]), indexBuffer);
            }

            for (int i = 0; i < 12; i++)
            {
                HoctreeNode[] faceNodes = new HoctreeNode[2];
                int[] c = { cellProcFaceMask[i][0], cellProcFaceMask[i][1] };

                faceNodes[0] = GetNode(childs[c[0]]);
                faceNodes[1] = GetNode(childs[c[1]]);

                ContourFaceProc(faceNodes, cellProcFaceMask[i][2], indexBuffer);
            }

            for (int i = 0; i < 6; i++)
            {
                HoctreeNode[] edgeNodes = new HoctreeNode[4];
                int[] c = new int[4]
                {
                    cellProcEdgeMask[i][0],
                    cellProcEdgeMask[i][1],
                    cellProcEdgeMask[i][2],
                    cellProcEdgeMask[i][3],
                };

                for (int j = 0; j < 4; j++)
                {
                    edgeNodes[j] = GetNode(childs[c[j]]);
                }

                ContourEdgeProc(edgeNodes, cellProcEdgeMask[i][4], indexBuffer);
            }
        }
    }

    #endregion

    #region Hashed Octree functions

    //Not much to say. I will probably remember this.. probably

    public short GetDepthFromNode(HoctreeNode node)
    {
        uint lc = node.LocCode;
        //11 is the max size for a uint(32bit)
        for(short i = 0; i < 11; ++i)
        {
            if (lc == 1)
                return i;
            lc = (lc >> 3);
        }

        return -1; //Borked, no depth value found
    }

    public int[] GetDepthAndSize(HoctreeNode node)
    {
        uint lc = node.LocCode;
        int depth = 0;
        int estSize = 1; //Estimated minimum size (always 1)

        //11 is the max size for a uint(32bit)
        for (int i = 0; i < 11; ++i)
        {
            if (lc == 1)
            {
                depth = i;
                //If size was not correct, do a reverse loop to figure out size
                if(estSize != Size)
                {
                    estSize = Size;
                    for(int j = 0; j < depth; ++j)
                    {
                        estSize /= 2;
                    }
                }
                return new int[] { depth, estSize };
            }
            else
            {
                lc = (lc >> 3);
                estSize *= 2;
            }
        }

        //broken, no clue
        return new int[] { -1, -1 };
    }

    public Vector3 GetMinFromNode(HoctreeNode node)
    {
        uint lc = node.LocCode;
        List<Vector3> mins = new List<Vector3>();

        //Calculate depth first, since we need it for size calc
        for (int i = 0; i < 11; ++i)
        {
            //When we have reached the top node, start reversing the node min list
            if (lc == 1)
            {
                //Parent will always be at 0,0,0 so no need to dig deeper
                Vector3 min = Position;
                int nSize = Size / 2;

                for (int j = mins.Count - 1; j >= 0; --j)
                {
                    min += mins[j] * nSize;
                    nSize /= 2;
                }
                return min;
            }
            else
            {
                mins.Add(CHILD_MIN_OFFSETS[(lc & 7)]);
                lc = (lc >> 3);
            }
        }
        return Vector3.zero; //Borked, no depth value found
    }


    public HoctreeNode[] GetChildsFromNode(HoctreeNode node)
    {
        HoctreeNode[] childs = new HoctreeNode[8];
        if (node.Type == NODE_TYPE.LEAF)
            return childs;

        for(uint i = 0; i < 8; ++i)
        {
            uint lc = node.LocCode << 3;
            lc += i;
            childs[i] = GetNode(lc);
        }
        return childs;
    }

    public uint[] GetChildLocsFromNode(HoctreeNode node)
    {
        uint[] childs = new uint[8];
        if (node.Type == NODE_TYPE.LEAF)
            return childs;

        for (uint i = 0; i < 8; ++i)
        {
            childs[i] = (node.LocCode << 3) + i;
        }
        return childs;
    }

    #endregion


    #region Static variables

    //Numbers.. A lot of this is from prior ported version of NGildea's work.
    //In short it's look up tables for which node should be used, vertices etc..

    public static int MATERIAL_AIR = 0;
    public static int MATERIAL_SOLID = 1;

    public static float QEF_ERROR = 1e-6f;
    public static int QEF_SWEEPS = 4;


    public static readonly Vector3[] CHILD_MIN_OFFSETS =
    {
	    // needs to match the vertMap from Dual Contouring impl
	    new Vector3( 0, 0, 0 ),
        new Vector3( 0, 0, 1 ),
        new Vector3( 0, 1, 0 ),
        new Vector3( 0, 1, 1 ),
        new Vector3( 1, 0, 0 ),
        new Vector3( 1, 0, 1 ),
        new Vector3( 1, 1, 0 ),
        new Vector3( 1, 1, 1 ),
    };

    public static readonly int[][] vertMap = new int[8][]
    {
        new int[3]{0,0,0},
        new int[3]{0,0,1},
        new int[3]{0,1,0},
        new int[3]{0,1,1},
        new int[3]{1,0,0},
        new int[3]{1,0,1},
        new int[3]{1,1,0},
        new int[3]{1,1,1}
    };

    // data from the original DC impl, drives the contouring process

    public static readonly int[][] edgevmap = new int[12][]
    {
        new int[2]{2,4},new int[2]{1,5},new int[2]{2,6},new int[2]{3,7},	// x-axis 
	    new int[2]{0,2},new int[2]{1,3},new int[2]{4,6},new int[2]{5,7},	// y-axis
	    new int[2]{0,1},new int[2]{2,3},new int[2]{4,5},new int[2]{6,7}		// z-axis
    };

    public static readonly int[] edgemask = { 5, 3, 6 };

    public static readonly int[][] faceMap = new int[6][]
    {
        new int[4]{4, 8, 5, 9},
        new int[4]{6, 10, 7, 11},
        new int[4]{0, 8, 1, 10},
        new int[4]{2, 9, 3, 11},
        new int[4]{0, 4, 2, 6},
        new int[4]{1, 5, 3, 7}
    };

    public static readonly int[][] cellProcFaceMask = new int[12][]
    {
        new int[3]{0,4,0},
        new int[3]{1,5,0},
        new int[3]{2,6,0},
        new int[3]{3,7,0},
        new int[3]{0,2,1},
        new int[3]{4,6,1},
        new int[3]{1,3,1},
        new int[3]{5,7,1},
        new int[3]{0,1,2},
        new int[3]{2,3,2},
        new int[3]{4,5,2},
        new int[3]{6,7,2}
    };

    public static readonly int[][] cellProcEdgeMask = new int[6][]
    {
        new int[5]{0,1,2,3,0},
        new int[5]{4,5,6,7,0},
        new int[5]{0,4,1,5,1},
        new int[5]{2,6,3,7,1},
        new int[5]{0,2,4,6,2},
        new int[5]{1,3,5,7,2}
    };

    public static readonly int[][][] faceProcFaceMask = new int[3][][]
    {
        new int[4][]{ new int[3]{4,0,0}, new int[3]{5,1,0}, new int[3]{6,2,0}, new int[3]{7,3,0} },
        new int[4][]{ new int[3]{2,0,1}, new int[3]{6,4,1}, new int[3]{3,1,1}, new int[3]{7,5,1} },
        new int[4][]{ new int[3]{1,0,2}, new int[3]{3,2,2}, new int[3]{5,4,2}, new int[3]{7,6,2} }
    };

    public static readonly int[][][] faceProcEdgeMask = new int[3][][]
    {
        new int[4][]{new int[6]{1,4,0,5,1,1},new int[6]{1,6,2,7,3,1},new int[6]{0,4,6,0,2,2},new int[6]{0,5,7,1,3,2}},
        new int[4][]{new int[6]{0,2,3,0,1,0},new int[6]{0,6,7,4,5,0},new int[6]{1,2,0,6,4,2},new int[6]{1,3,1,7,5,2}},
        new int[4][]{new int[6]{1,1,0,3,2,0},new int[6]{1,5,4,7,6,0},new int[6]{0,1,5,0,4,1},new int[6]{0,3,7,2,6,1}}
    };

    public static readonly int[][][] edgeProcEdgeMask = new int[3][][]
    {
        new int[2][]{new int[5]{3,2,1,0,0},new int[5]{7,6,5,4,0}},
        new int[2][]{new int[5]{5,1,4,0,1},new int[5]{7,3,6,2,1}},
        new int[2][]{new int[5]{6,4,2,0,2},new int[5]{7,5,3,1,2}},
    };

    public static readonly int[][] processEdgeMask = new int[3][]
    {
        new int[4]{3,2,1,0},new int[4]{7,5,6,4},new int[4]{11,10,9,8}
    };

    #endregion
}
