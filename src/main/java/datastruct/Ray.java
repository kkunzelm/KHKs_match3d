package datastruct;/*
					Copyright 2008, 2009, 2010  Peter Hofmann
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
 * Ein einzelner datastruct.Ray.
 */
public class Ray {

	private static final double corrEps = 1e-7;

	public final Vec3 origin;
	public final Vec3 direction;

	// KHK nicht nötig: public datastruct.Vec3 reciDir;

	/**
	 * Erzeugt einen neuen datastruct.Ray an dieser Position mit dieser Richtung
	 */
	public Ray(Vec3 origin, Vec3 direction) {
		this.origin = new Vec3(origin);
		this.direction = direction;
		this.direction.normalize();

		// Um ein kleines Epsilon schubsen - das ist nötig, um
		// Rundungsfehler auszugleichen
		this.origin.x += (this.direction.x * corrEps);
		this.origin.y += (this.direction.y * corrEps);
		this.origin.z += (this.direction.z * corrEps);

		// Kehrwerte für schnelleren AABB-datastruct.Ray-Test
		// Keine Angst vor Division durch 0: Dann kommt +/- infinity raus,
		// was auch dort im Test für den Größenvergleich benötigt wird.
		// KHK nicht nötig: reciDir = new datastruct.Vec3(1.0 / direction.x, 1.0 /
		// direction.y, 1.0 / direction.z);
	}

	/**
	 * Gibt den konkreten Ort zurück, der sich "alpha" weit vom datastruct.Ray-
	 * Ursprung entfernt befindet.
	 */
	public Vec3 evaluate(double alpha) {
		Vec3 out = new Vec3(origin);
		out.x += (direction.x * alpha);
		out.y += (direction.y * alpha);
		out.z += (direction.z * alpha);
		return out;
	}

	public String toString() {
		return "datastruct.Ray(new " + origin + ", new " + direction + ")";
	}
}