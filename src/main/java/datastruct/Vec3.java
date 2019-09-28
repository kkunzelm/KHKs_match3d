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

import java.util.Random;

/**
 * Selbsterklärende Klasse mit Vektor-Operationen
 */
public class Vec3
{
	public double x, y, z;

	public Vec3(double x, double y, double z)
	{
		this.x = x;
		this.y = y;
		this.z = z;
	}

	public Vec3(Vec3 other)
	{
		this(other.x, other.y, other.z);
	}

	public Vec3()
	{
		this(0, 0, 0);
	}

	/**
	 * Nur damit via i über die Achsen iteriert werden kann.
	 */
	public double getAxis(int which)
	{
		switch (which)
		{
			case 0:
				return x;
			case 1:
				return y;
			case 2:
				return z;
			default:
				return 0.0;
		}
	}

	/**
	 * Nur damit via i über die Achsen iteriert werden kann.
	 */
	public void setAxis(int which, double val)
	{
		switch (which)
		{
			case 0:
				x = val;
				break;
			case 1:
				y = val;
				break;
			case 2:
				z = val;
				break;
		}
	}

	public void add(Vec3 other)
	{
		x += other.x;
		y += other.y;
		z += other.z;
	}

	public void subtract(Vec3 other)
	{
		x -= other.x;
		y -= other.y;
		z -= other.z;
	}

	public Vec3 plus(Vec3 other)
	{
		return new Vec3(x + other.x, y + other.y, z + other.z);
	}

	public Vec3 minus(Vec3 other)
	{
		return new Vec3(x - other.x, y - other.y, z - other.z);
	}

	public void scale(double a)
	{
		x *= a;
		y *= a;
		z *= a;
	}

	public Vec3 times(double a)
	{
		Vec3 out = new Vec3(this);
		out.scale(a);
		return out;
	}

	public double dot(Vec3 other)
	{
		return x * other.x + y * other.y + z * other.z;
	}

	public Vec3 cross(Vec3 other)
	{
		return new Vec3(
							y * other.z - z * other.y,
							z * other.x - x * other.z,
							x * other.y - y * other.x
						);
	}

	public double length()
	{
		return Math.sqrt(dot(this));
	}

	public double lengthSquared()
	{
		return dot(this);
	}

	public void normalize()
	{
		scale(1.0 / length());
	}

	public Vec3 normalized()
	{
		return times(1.0 / length());
	}

	public Vec3 jittered(double epsilon, Random rGen)
	{
		// nextGaussian() sieht besser aus, ist aber definitiv langsamer.
		return new Vec3(this.x + (rGen.nextGaussian() * 0.5) * epsilon,
						this.y + (rGen.nextGaussian() * 0.5) * epsilon,
						this.z + (rGen.nextGaussian() * 0.5) * epsilon
						);
		/*
		return new datastruct.Vec3(this.x + (Math.random() - 0.5) * epsilon,
						this.y + (Math.random() - 0.5) * epsilon,
						this.z + (Math.random() - 0.5) * epsilon
						);
		*/
	}

	public String toString()
	{
		return "datastruct.Vec3(" + x + ", " + y + ", " + z + ")";
	}
}