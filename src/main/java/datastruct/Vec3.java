package datastruct;// KHK: status 30.4.2014
// KHK: I downloaded this file with the accompaning files 
// KHK: I used it at an intermediate devolopment state for my target to source transformation
//      it worked fine, but was very slow
//      I learned a lot from these files. Especially the math behind raytracing!! Many thanks to the author.
//      I also like the data structures he implemented: vec3, Triangle3d...
//      This implementation of Triangle3d was modified by me to add several methods to return values (like boundingBox ...)

/*
	Copyright 2008, 2009, 2010  Peter Hofmann
        http://www.uninformativ.de

        http://www.uninformativ.de/?section=news&ndo=single&newsid=117
        http://www.uninformativ.de/bin/RaytracingSchnitttests-76a577a-CC-BY.pdf

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful, but
	WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
	General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
 * Selbsterklärende Klasse mit Vektor-Operationen
 */
public class Vec3 {
	public double x, y, z;

	public Vec3(double x, double y, double z) {
		this.x = x;
		this.y = y;
		this.z = z;
	}

	public Vec3(Vec3 other) {
		this(other.x, other.y, other.z);
	}

	public Vec3() {
		this(0, 0, 0);
	}

	public Vec3 minus(Vec3 other) {
		return new Vec3(x - other.x, y - other.y, z - other.z);
	}

	private void scale(double a) {
		x *= a;
		y *= a;
		z *= a;
	}

	public Vec3 times(double a) {
		Vec3 out = new Vec3(this);
		out.scale(a);
		return out;
	}

	public double dot(Vec3 other) {
		return x * other.x + y * other.y + z * other.z;
	}

	public Vec3 cross(Vec3 other) {
		return new Vec3(y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x);
	}

	private double length() {
		return Math.sqrt(dot(this));
	}

	public void normalize() {
		scale(1.0 / length());
	}

	public String toString() {
		return "datastruct.Vec3(" + x + ", " + y + ", " + z + ")";
	}
}