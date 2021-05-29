using System.Collections.Generic;
using UnityEngine;

public struct MeshVertex
{
	public Vector3 xyz, normal;
	public Color color;
	public MeshVertex(Vector3 _xyz, Vector3 _normal, Color _color)
    {
		xyz = _xyz;
		normal = _normal;
		color = _color;
	}
}
