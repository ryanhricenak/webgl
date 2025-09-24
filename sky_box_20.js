"use strict";

var canvas;
var gl;

var numVertices = 36;
var deg_to_rad = Math.PI / 180;

var pointsArray = [];
var colorsArray = [];

// sizes of different objects
var sky_start;
var mountain_start;
var water_start;
var grass_start;
var sand_start;
var tree1_start;
var tree2_start;
var tree3_start;
var tree4_start;
var tree5_start;
var tree6_start;

var pi = Math.PI;

var wf = 200;

var water_trig = [];

var vertices = [
    vec4(-0.5, -0.5, 0.5, 1.0),
    vec4(-0.5, 0.5, 0.5, 1.0),
    vec4(0.5, 0.5, 0.5, 1.0),
    vec4(0.5, -0.5, 0.5, 1.0),
    vec4(-0.5, -0.5, -0.5, 1.0),
    vec4(-0.5, 0.5, -0.5, 1.0),
    vec4(0.5, 0.5, -0.5, 1.0),
    vec4(0.5, -0.5, -0.5, 1.0),
];

var vertexColors = [
    vec4(0.0, 0.0, 0.0, 1.0), // black
    vec4(1.0, 0.0, 0.0, 1.0), // red
    vec4(1.0, 1.0, 0.0, 1.0), // yellow
    vec4(0.0, 1.0, 0.0, 1.0), // green
    vec4(0.0, 0.0, 1.0, 1.0), // blue
    vec4(1.0, 0.0, 1.0, 1.0), // magenta
    vec4(0.0, 1.0, 1.0, 1.0), // cyan
    vec4(1.0, 1.0, 1.0, 1.0), // white
];

var near = -0.03;
var far = 0.03;
var radius = 0.56;
var theta = 0.35; //45.0 * Math.PI/180.0;
var phi = 270.0 * Math.PI / 180.0;
var dr = 1 * Math.PI / 180.0;

var left = -1.0;
var right = 1.0;
var ytop = 1.0;
var bottom = -1.0;

// save terrian data points
var data = [];
var TerrainHi = 0.25
var terrRows = 512;
var terrCols = 512;
var index = 0;

// save into 1D Array

var mvMatrix, pMatrix;
var modelView, projection;
var eye;

var at = vec3(0.0, 0.0, 0.0); /*look at origin*/
const up = vec3(0.0, 1.0, 0.0); /*y-axis direction as up*/

var pause_flag = 1;

// code for mountain

var mountain_center = vec3(-0.5, 0.0, -0.5); // center of the base of the mountain
var mountain_height = 3 / 4;
var mountain_radius = 0.5;

var light = vec3(0, 3.0, 0); // light source for the mountain

var mountain_red = 134 / 255;
var mountain_blue = 136 / 255;
var mountain_green = 139 / 255;

var mountain_normal_vectors = [];

var sides = 40;
var num_rows_per_triangle = 40;
var points_per_triangle;

var water_normal_vectors = [];
//var water_trigs = [];

// Generate terrian data using trig functions
function genTerrainData(rows, cols) {
    var xinc = 2 / rows;
    var zinc = 2 / cols;

    var x = -1;

    for (let i = 0; i < rows; x += xinc, i++) {
        data[i] = new Array(cols);
        var z = -1;
        for (let j = 0; j < cols; z += zinc, j++) {
            if (j <= (9 * cols / 20)) {
                data[i][j] = Math.abs(((Math.cos(3 * j) + Math.sin(10 * i)) / 20)) + (TerrainHi / 4); // this is the grasslands for the terrain			
            } else if ((9 * cols / 20) <= j && j <= 3 * cols / 5) {
                data[i][j] = -0.5; //Math.cos(j * Math.PI / 120) / 20 + (TerrainHi / 4); // this is the sand for the terrain
            } else {
                data[i][j] = -0.5;  // Math.cos(j * Math.PI / 120) / 20 + (TerrainHi / 4) // this is the water 
            }
        }
    }
}

// Prepare mesh data to bufferData
function prepMesh(nRows, nColumns) {

    for (var i = 0; i < nRows; ++i) {
        for (var j = 0; j < nColumns; ++j) {
            pointsArray[index] = vec4(2 * i / nRows - 1, data[i][j], 2 * j / nColumns - 1, 1.0);
            //colorsArray[index] = vec4(i/nRows, j/nColumns, 0.0, 1.0);

            if (j <= (9 * nColumns / 20)) {
                colorsArray[index] = vec4(0.0, j / nColumns, 0.0, 1.0); // green for grass
            } else if (((9 * nColumns / 20) <= j) && (j <= (3 * nColumns / 5))) {
                colorsArray[index] = vec4(j / nRows, j / nColumns, 0.0, 1.0); // yellow for sand
            } else {
                colorsArray[index] = vec4(0.0, 0.0, j / nColumns, 1.0); // blue for water
            }

            index++;
        }
    }

    for (var j = 0; j < nColumns; ++j) {
        for (var i = 0; i < nRows; ++i) {
            pointsArray[index] = vec4(2 * i / nRows - 1, data[i][j], 2 * j / nColumns - 1, 1.0);
            colorsArray[index] = vec4(0.0, 0.0, 0.0, 1.0);
            index++;
        }
    }

    sky_start = pointsArray.length;
}

function generate_mountain() {

    var side_inc = 360 / sides;
    var flag = 1;

    for (let i = 0; i <= side_inc * sides; i += side_inc) {

        var p1 = vec4(mountain_center[0], mountain_height, mountain_center[2], 1.0);
        var p2 = vec4(mountain_center[0] + mountain_radius * Math.cos(i * deg_to_rad),				mountain_center[1], mountain_center[2] + mountain_radius * Math.sin(i * deg_to_rad), 1.0);
        var p3 = vec4(mountain_center[0] + mountain_radius * Math.cos((i + side_inc) * deg_to_rad), mountain_center[1], mountain_center[2] + mountain_radius * Math.sin((i + side_inc) * deg_to_rad), 1.0);

        generate_unit_triangle(p1, p2, p3, num_rows_per_triangle);

        if (flag > 0) {
            points_per_triangle = pointsArray.length - mountain_start;
            flag *= (-1);
        }
    }
	apply_colors();
}

function generate_unit_triangle(a, b, c, d) {
    var p1 = vec4(a[0] / d, a[1] / d, a[2] / d, 1.0)
    var p2 = vec4(b[0] / d, b[1] / d, b[2] / d, 1.0)
    var p3 = vec4(c[0] / d, c[1] / d, c[2] / d, 1.0)

    let q = b[0] - p2[0];
    let w = b[1] - p2[1];
    let r = b[2] - p2[2];

    var diff = vec4(q, w, r, 1.0);

    generate_triangle(p1, p2, p3, d, diff);
}

var normal_vector_flag = 1;
var first_triangle_in_row_flag = 1;

var final_c_in_row = [];
var final_c_counter = 0;

var first_triangle_flag = 1;

var b_array = [];
var b_count = 0;

var original_a;
var original_b;
var original_c;

function generate_triangle(a, b, c, d, diff) {
    var cb_x_inc = c[0] - b[0];
    var cb_y_inc = c[1] - b[1];
    var cb_z_inc = c[2] - b[2];

    var ab_x_inc = a[0] - b[0];
    var ab_y_inc = a[1] - b[1];
    var ab_z_inc = a[2] - b[2];

    // start with the bottom
    let e = d;

    let q = diff[0];
    let w = diff[1];
    let r = diff[2];

    a = vec4(a[0] + q, a[1] + w, a[2] + r, 1.0);
    b = vec4(b[0] + q, b[1] + w, b[2] + r, 1.0);
    c = vec4(c[0] + q, c[1] + w, c[2] + r, 1.0);

	var last_c =  vec4( c[0] + (cb_x_inc * num_rows_per_triangle), c[1] + (cb_y_inc * num_rows_per_triangle), c[2] + (cb_z_inc * num_rows_per_triangle), 1.0 );

    var normal_a = a;
    var normal_b = b;
    var normal_c = c;

    var u;
    var v;

    var old_a_array = [];
    var old_a_count = 0;
    var first_row_flag = 1;

	if(first_triangle_flag > 0){
		// poop;
	} else{
		b = original_c;
	}
	
	original_c = c;
	original_a = a;

    for (let j = 0; j < d; ++j) {

		if(first_row_flag > 0){
			pointsArray.push(vec4(b[0], b[1], b[2], 1.0));		
		} else{
			b = old_a_array[old_a_count];
			++old_a_count;
			pointsArray.push(b);
		}
		
        let cx = 0;
        let cz = 0;

        u = vec4(b[0] - a[0], b[1] - a[1], b[2] - a[2], 1.0);
        v = vec4(c[0] - b[0], c[1] - b[1], c[2] - b[2], 1.0);

        calc_normal_vector(u, v); // three because we need all the points from the first triangle 
        calc_normal_vector(u, v);
        calc_normal_vector(u, v);

        var old_a;
        var old_c;

        for (let i = 0; i < e; ++i) {
			
			if(i == 0){
				normal_a = original_a;
			} else{
				normal_a = vec4(a[0] + cx + mountain_bump_random(), a[1], a[2] + cz + mountain_bump_random(), 1.0);				
			}

			if ( (i + 1) == e ){
				normal_c = last_c;
				++old_a_count;
			} else if (first_row_flag < 0) {
                normal_c = old_a_array[old_a_count];
                ++old_a_count;
            } else {
                normal_c = vec4(c[0] + cx + mountain_bump_random(), c[1], c[2] + cz + mountain_bump_random(), 1.0);
            }

            pointsArray.push(normal_a);
            pointsArray.push(normal_c);

            if (first_triangle_in_row_flag < 0) {
                u = vec4(old_a[0] - normal_a[0], old_a[1] - normal_a[1], old_a[2] - normal_a[2], 1.0);
                v = vec4(old_c[0] - old_a[0], old_c[1] - old_a[1], old_c[2] - old_a[2], 1.0);

                calc_normal_vector(u, v);

                u = vec4(old_c[0] - normal_a[0], old_c[1] - normal_a[1], old_c[2] - normal_a[2], 1.0);
                v = vec4(normal_c[0] - old_c[0], normal_c[1] - old_c[1], normal_c[2] - old_c[2], 1.0);

                calc_normal_vector(u, v);
            } else {
                first_triangle_in_row_flag *= (-1);
            }

            cx += cb_x_inc;
            cz += cb_z_inc;

            old_a_array.push(normal_a);

            old_a = normal_a;
            old_c = normal_c;
        }

		b = vec4(a[0], a[1], a[2], 1.0);
        c = vec4(b[0] + cb_x_inc, b[1] + cb_y_inc, b[2] + cb_z_inc, 1.0);
        a = vec4(b[0] + ab_x_inc, b[1] + ab_y_inc, b[2] + ab_z_inc, 1.0);

        --e;
        first_triangle_in_row_flag *= (-1);
        first_row_flag = -1;
    }
	
    first_triangle_flag = -3;
}

function calc_normal_vector(u, v) {
    mountain_normal_vectors.push(vec3((u[1] * v[2]) - (u[2] * v[1]), (u[2] * v[0]) - (u[0] * v[2]),
        (u[0] * v[1]) - (u[1] * v[0])));
}

function apply_colors() {
	
    for (let o = 0; o < mountain_normal_vectors.length; ++o) {
        let i = o + mountain_start;

        var n = mountain_normal_vectors[o];
        var p = pointsArray[i];
		
        var l = vec4(n[0] - p[0], n[1] - p[1], n[2] - p[2], 1.0);

        var n_mag = Math.sqrt((n[0] ** 2) + (n[1] ** 2) + (n[2] ** 2));
        var l_mag = Math.sqrt((l[0] ** 2) + (l[1] ** 2) + (l[2] ** 2));

        var beta = ((n[0] * l[0]) + (n[1] * l[1]) + (n[2] * l[2])) / (n_mag * l_mag);

        if (beta < 0) {
            beta = 0;
        }

        colorsArray.push(vec4(mountain_red * beta, mountain_green * beta, mountain_blue * beta, 1.0));
    }
}

function generate_grass_underground(){
	let a = 1 / 10;
	var p1 = vec4( 1.0, a, -1.0, 1.0);
	var p2 = vec4(-1.0, a, -1.0, 1.0);
	var p3 = vec4( 1.0, a, -0.1, 1.0);
	var p4 = vec4(-1.0, a, -0.1, 1.0);

	pointsArray.push(p1);
	pointsArray.push(p2);
	pointsArray.push(p3);
	pointsArray.push(p4);	
	
	colorsArray.push( vec4(0.0, 0.3, 0.0, 1.0) );
	colorsArray.push( vec4(0.0, 0.3, 0.0, 1.0) );
	colorsArray.push( vec4(0.0, 0.3, 0.0, 1.0) );
	colorsArray.push( vec4(0.0, 0.3, 0.0, 1.0) );
}

function water_rand(){
	let min = -1;
	let max =  1;
	return ( Math.random() * (max - min + 1) + min ) / 1000;
}

function random_deg(){
	let min = 0;
	let max = 90;
	return (Math.random() * (max - min + 1) + min);
}

function generate_trig(){
	for(let i = 0; i < wf; ++i){
		water_trig.push( ( (wf - i) * pi) / wf );
	}
	
	for(let i = 0; i < wf; ++i){
		water_trig.push( i * pi / wf );
	}
}

var water_row;
var water_row_flag = 1;
var water_length;

var rows_of_water = 90;
var cols_of_water = 90;

var top_water = [];
var bot_water = [];

var first_water_row_flag = 1;
var first_water_normal_vector_flag = 1;

var trig_flag = 1;
var trig_index = 0;

function generate_water_first_time(){
	let h = 1 / 20;

	var old_a;
	var old_b;

	water_length = 0;

	var i_inc = 2 / rows_of_water;
	var j_inc = (0.8) / cols_of_water;

	var a;
	var b;
	var u;
	var v;

	for(let i = 0; i < rows_of_water; ++i){
		for(let j = 0; j < cols_of_water; ++j){
			if(first_water_row_flag > 0){
				a = vec4( -1.0 + ( i * i_inc ), h + water_rand(), j * j_inc + 0.2, 1.0);
				b = vec4( -1.0 + ( (i + 1) * i_inc), h + water_rand(), j * j_inc + 0.2, 1.0);			
								
				if(j === 0){
					// basically just chill for a moment
				} else{
					// calculate the normal vectors for a and b
					u = vec4( Math.abs(old_a[0] - a[0]), Math.abs(old_a[1] - a[1]), Math.abs(old_a[2] - a[2]), 1.0 );
					v = vec4( Math.abs(old_b[0] - old_a[0]), Math.abs(old_b[1] - old_a[1]), Math.abs(old_b[2] - old_a[2]), 1.0 );
					
					if(first_water_normal_vector_flag > 0){
						calc_normal_vector_water(u, v);
						calc_normal_vector_water(u, v);
						calc_normal_vector_water(u, v);
						first_water_normal_vector_flag = -1;
					} else{
						calc_normal_vector_water(u, v);
					}
					
					u = vec4( Math.abs(old_b[0] - a[0]), Math.abs(old_b[1] - a[1]), Math.abs(old_b[2] - a[2]), 1.0 );
					v = vec4( Math.abs(b[0] - old_b[0]), Math.abs(b[1] - old_b[1]), Math.abs(b[2] - old_b[2]), 1.0 );
					calc_normal_vector_water(u, v);
					}
				
				bot_water.push(b);
				
				pointsArray.push(a);
				pointsArray.push(b);
				water_length += 2;
				
				old_a = a;
				old_b = b;
				
			} else{
				a = top_water[j];
				b = vec4( -1.0 + ( (i + 1) * i_inc) , h + water_rand(), j * j_inc + 0.2, 1.0);
							
				if(j === 0){
					// chillax
				} else{
										// calculate the normal vectors for a and b
					u = vec4( Math.abs(old_a[0] - a[0]), Math.abs(old_a[1] - a[1]), Math.abs(old_a[2] - a[2]), 1.0 );
					v = vec4( Math.abs(old_b[0] - old_a[0]), Math.abs(old_b[1] - old_a[1]), Math.abs(old_b[2] - old_a[2]), 1.0 );
					
					if(first_water_normal_vector_flag > 0){
						calc_normal_vector_water(u, v);
						calc_normal_vector_water(u, v);
						calc_normal_vector_water(u, v);
						first_water_normal_vector_flag = -1;
					} else{
						calc_normal_vector_water(u, v);
					}
										
					u = vec4( Math.abs(old_b[0] - a[0]), Math.abs(old_b[1] - a[1]), Math.abs(old_b[2] - a[2]), 1.0 );
					v = vec4( Math.abs(b[0] - old_b[0]), Math.abs(b[1] - old_b[1]), Math.abs(b[2] - old_b[2]), 1.0 );
					calc_normal_vector_water(u, v);
					
				}
				
				pointsArray.push(a);
				pointsArray.push(b);
				water_length += 2;
				bot_water.push(b);
			}
		}

		first_water_normal_vector_flag = 1;
		first_water_row_flag = -1;
		top_water = bot_water;
		bot_water = [];
	}

	apply_water_colors();
}

var water_trig_index = 0;

var check_flaf = 1;

function generate_water(){
	var points = [];
	// go through the list and get the new heights
		
	let w = water_start;
	var index;
	
	var p;
	index = 0;
	
	for(let i = 0; i < 2 * rows_of_water; ++i){
		for(let j = 0; j < cols_of_water; ++j){
			index = ( (i * rows_of_water) + j ) + water_start;
			
			p = pointsArray[index];
			
			
			let a = vec4(p[0], p[1] - Math.sin(water_trig[ ( (cols_of_water-j) + water_trig_index) % wf ]) / 100, p[2], 1.0 );
			points.push( a );

		}
	}

	++water_trig_index;
	if(water_trig_index === wf){
		water_trig_index = 0;
	}
	
	check_flaf = -1;
	
	// ++water_trig_index;
	
	// now that we have all the points, sub buffer my guy
	gl.bindBuffer(gl.ARRAY_BUFFER, vBuffer);
    gl.bufferSubData(gl.ARRAY_BUFFER, (water_start * 16), flatten(points));
}

function calc_normal_vector_water(u, v){
		let a = (u[1] * v[2]) - (u[2] * v[1]);
		let b = (u[2] * v[0]) - (u[0] * v[2]);
		let c = (u[0] * v[1]) - (u[1] * v[0]);
		
	    water_normal_vectors.push(vec3( a, b, c) );
}

var check_flag = 1;

function apply_water_colors(){
	var water_red = 15 / 255;
	var water_green = 94 / 255;
	var water_blue = 156 / 255;

	for (let o = 0; o < water_normal_vectors.length; ++o) {
        let i = o + water_start;

        var n = water_normal_vectors[o];
		
		var p = pointsArray[i];
        var l = vec4(light[0] - p[0], light[1] - p[1], light[2] - p[2], 1.0);
	
        var n_mag = Math.sqrt((n[0] ** 2) + (n[1] ** 2) + (n[2] ** 2));
        var l_mag = Math.sqrt((l[0] ** 2) + (l[1] ** 2) + (l[2] ** 2));

        var beta = ((n[0] * l[0]) + (n[1] * l[1]) + (n[2] * l[2])) / (n_mag * l_mag);
		
		beta = Math.abs(beta);
		if( Number.isNaN(beta) ){
			beta = 0;
		}
		
        colorsArray.push(vec4(water_red * beta, water_green * beta, water_blue * beta, 1.0));
    }
}

function mountain_bump_random() {
    let min = -5;
    let max = 5;
    return (Math.random() * (max - min + 1) + min) / 1000;
}

var sand_length;
var rows_of_sand = 20;
var cols_of_sand = 20;

var first_sand_row_flag = 1;
var first_sand_normal_vector_flag = 1;

var sand_normal_vectors = [];

function calc_normal_vector_sand(u, v){
	    sand_normal_vectors.push(vec3((u[1] * v[2]) - (u[2] * v[1]), (u[2] * v[0]) - (u[0] * v[2]),
        (u[0] * v[1]) - (u[1] * v[0])));

}

function sand_rand(){
	var max = 1;
	var min = -1;
	return (Math.random() * (max - min + 1) + min) / 1000;
}

function generate_sand(){
	var sand = [];
	
	let h = 1 / 10;

	var old_a;
	var old_b;

	sand_length = 0;

	var i_inc = 2 / rows_of_sand;
	var j_inc = -(0.45) / cols_of_sand;

	var top_sand = [];
	var bot_sand = [];

	var a;
	var b;
	var u;
	var v;

	var p1 = vec4( -1.0, 0.10, -0.2, 1.0 );
	var p2 = vec4( -1.0, 0.0, 0.6, 1.0 ); // need p2 and 3 

	var p3 = vec4( 1.0, 0.10, -0.2, 1.0 );
	var p4 = vec4( 1.0, 0.0, 0.6, 1.0 );

	var do_it_once = 1;
	var do_index = 0;

	for(let i = 0; i < rows_of_sand; ++i){
		for(let j = 0; j < cols_of_sand; ++j){
			if(first_sand_row_flag > 0){
				a = vec4( -1.0 + ( i * i_inc ), h - ( (cols_of_sand - j) / 300) + sand_rand(), j * j_inc + 0.3, 1.0);
				b = vec4( -1.0 + ( (i + 1) * i_inc), h - ( (cols_of_sand - j) / 300) + sand_rand(), j * j_inc + 0.3, 1.0);			
								
				if(j === 0){
					// basically just chill for a moment
				} else{
					// calculate the normal vectors for a and b
					u = vec4( (old_a[0] - a[0]), old_a[1] - a[1], old_a[2] - a[2], 1.0 );
					v = vec4( (old_b[0] - old_a[0]), (old_b[1] - old_a[1]), (old_b[2] - old_a[2]), 1.0 );					
					
					if(first_sand_normal_vector_flag > 0){
						calc_normal_vector_sand(u, v);
						calc_normal_vector_sand(u, v);
						calc_normal_vector_sand(u, v);
						first_sand_normal_vector_flag = -1;
					} else{
						calc_normal_vector_sand(u, v);
					}
					
					u = vec4( (old_b[0] - a[0]), (old_b[1] - a[1]), (old_b[2] - a[2]), 1.0 );
					v = vec4( (b[0] - old_b[0]), (b[1] - old_b[1]), (b[2] - old_b[2]), 1.0 );
					calc_normal_vector_sand(u, v);
				}
				
				bot_sand.push(b);
				
				pointsArray.push(a);
				pointsArray.push(b);
				sand.push(1);
				sand.push(1);
				
				old_a = a;
				old_b = b;
				
			} else{
				a = top_sand[j];
				b = vec4( -1.0 + ( (i + 1) * i_inc) , h - ( (cols_of_sand - j) / 300) + sand_rand(), j * j_inc + 0.3, 1.0);
							
				if(j === 0){
					// chillax
				} else{
										// calculate the normal vectors for a and b
					u = vec4( (old_a[0] - a[0]), (old_a[1] - a[1]), (old_a[2] - a[2]), 1.0 );
					v = vec4( (old_b[0] - old_a[0]), (old_b[1] - old_a[1]), (old_b[2] - old_a[2]), 1.0 );
					
					if(first_sand_normal_vector_flag > 0){
						calc_normal_vector_sand(u, v);
						calc_normal_vector_sand(u, v);
						calc_normal_vector_sand(u, v);
						first_sand_normal_vector_flag = -1;
					} else{
						calc_normal_vector_sand(u, v);
					}
					
					u = vec4( (old_b[0] - a[0]), (old_b[1] - a[1]), (old_b[2] - a[2]), 1.0 );
					v = vec4( (b[0] - old_b[0]), (b[1] - old_b[1]), (b[2] - old_b[2]), 1.0 );
					calc_normal_vector_sand(u, v);
				}
				
				pointsArray.push(a);
				pointsArray.push(b);
				sand.push(1);
				sand.push(1);
				bot_sand.push(b);
			}
		}
		
		do_it_once = -1;
		
		first_sand_normal_vector_flag = 1;
		first_sand_row_flag = -1;
		top_sand = bot_sand;
		bot_sand = [];
	}

	apply_sand_colors();
	sand_length = sand.length;
}

function apply_sand_colors(){
	var sand_red = 226 / 255;
	var sand_green = 202 / 255;
	var sand_blue = 118 / 255; 
	
	
	for (let o = 0; o < sand_normal_vectors.length; ++o) {
        let i = o + sand_start;

        var n = sand_normal_vectors[o];
		
		var p = pointsArray[i];
        var l = vec4(- light[0] + p[0], - light[1] + p[1], - light[2] + p[2], 1.0);

        var n_mag = Math.sqrt((n[0] ** 2) + (n[1] ** 2) + (n[2] ** 2));
        var l_mag = Math.sqrt((l[0] ** 2) + (l[1] ** 2) + (l[2] ** 2));

		if(n_mag <= 0){
			n_mag = 1;
		}
		
		if(l_mag <= 0){
			l_mag = 1;
		}

        var beta = ((n[0] * l[0]) + (n[1] * l[1]) + (n[2] * l[2])) / (n_mag * l_mag);
		
		beta = Math.max(beta, 0); 
		
        colorsArray.push(vec4(sand_red * beta, sand_green * beta, sand_blue * beta, 1.0));
    }
}


function colorCube() {
    quad(1, 0, 3, 2);
    quad(2, 3, 7, 6);
    quad(3, 0, 4, 7);
    quad(6, 5, 1, 2);
    quad(4, 5, 6, 7);
    quad(5, 4, 0, 1);
}

function quad(a, b, c, d) {
    var vertices = [
        vec4(-1.0, -1.0, 1.0, 1.0),
        vec4(-1.0, 1.0, 1.0, 1.0),
        vec4(1.0, 1.0, 1.0, 1.0),
        vec4(1.0, -1.0, 1.0, 1.0),
        vec4(-1.0, -1.0, -1.0, 1.0),
        vec4(-1.0, 1.0, -1.0, 1.0),
        vec4(1.0, 1.0, -1.0, 1.0),
        vec4(1.0, -1.0, -1.0, 1.0)
    ];

    // We need to parition the quad into two triangles in order for
    // WebGL to be able to render it.  In this case, we create two
    // triangles from the quad indices

    //vertex color assigned by the index of the vertex

    var indices = [a, b, c, a, c, d];

    for (var i = 0; i < indices.length; ++i) {
        pointsArray.push(vertices[indices[i]]);
        colorsArray.push(vec4(2 / 255, 204 / 255, 254 / 255, 1.0));
    }
}

const num_tree_faces_root = 100;
const sides_per_tree = 30;
const num_rows_per_tree_triangle = 40;

var total_number_of_trees = 6;

var tree_center1 = vec3(0.5, 0.0, -0.5); // center of the base of the tree
var tree_height1 = random_tree_height(); //random_tree_height();
var tree_radius1 = random_tree_radius(tree_height1); //random_tree_radius(tree_height1);
var tree_root_height1 = tree_height1 / 3;

var tree_center2 = vec3(0.8, 0.0, -0.36); // center of the base of the tree
var tree_height2 = random_tree_height();
var tree_radius2 = random_tree_radius(tree_height2);
var tree_root_height2 = tree_height2 / 3;

var tree_center3 = vec3(0.32, 0.0, -0.2); // center of the base of the tree
var tree_height3 = random_tree_height();
var tree_radius3 = random_tree_radius(tree_height3);
var tree_root_height3 = tree_height3 / 3;

var tree_center4 = vec3(0.72, 0.0, -0.89); // center of the base of the tree
var tree_height4 = random_tree_height();
var tree_radius4 = random_tree_radius(tree_height4);
var tree_root_height4 = tree_height4 / 3;

var tree_center5 = vec3(0.13, 0.0, -0.287); // center of the base of the tree
var tree_height5 = random_tree_height();
var tree_radius5 = random_tree_radius(tree_height5);;
var tree_root_height5 = tree_height5 / 3;

var tree_center6 = vec3(0.368, 0.0, -0.78); // center of the base of the tree
var tree_height6 = random_tree_height;
var tree_radius6 = random_tree_radius(tree_height6);;
var tree_root_height6 = tree_height6 / 3;

function random_tree_height(){
	var max = 2;
	var min = 1;
	return ( ( Math.random() * (max - min + 1) ) + min ) / 10;
}

function random_tree_radius(t){
	var max = 1.5 * t;
	var min = t;
	return ( ( Math.random() * (max - min + 1) ) + min ) / 10;
}

var points_per_tree_triangle;

var normal_vector_flag_tree = 1;
var first_triangle_in_tree_row_flag = 1;
var tree_normal_vectors = [];

var tree_red = 42 / 255;
var tree_green = 126 / 255;
var tree_blue = 25 / 255;

var tree_root_red = 210 / 255;
var tree_root_green = 113 / 255;
var tree_root_blue = 64 / 255;

function generate_tree(tree_center, tree_height, tree_radius, tree_root_height){
	
	var side_inc = 360 / sides_per_tree;
	var flag = 1;

	generate_tree_root(tree_center, tree_height, tree_radius, tree_root_height);

	for(let i = 0; i <= side_inc * sides_per_tree; i += side_inc){
	
		var p1 = vec4( tree_center[0], tree_height + tree_center[1] + tree_root_height, tree_center[2], 1.0);	
	
		var p2 = vec4(tree_center[0] + tree_radius * Math.cos(i * deg_to_rad), tree_center[1] + tree_root_height, tree_center[2] + tree_radius * Math.sin(i * deg_to_rad), 1.0);
		var p3 = vec4(tree_center[0] + tree_radius * Math.cos((i + side_inc) * deg_to_rad), tree_center[1] + tree_root_height, tree_center[2] + tree_radius * Math.sin((i + side_inc) * deg_to_rad), 1.0);

		var color_code = vec4(Math.random(), Math.random(), Math.random(), 1.0);
	
		generate_unit_triangle_tree(p1, p2, p3, num_rows_per_tree_triangle);
		
		if(flag > 0){
			points_per_tree_triangle = pointsArray.length - tree1_start;
			flag *= (-1);
		}		
	}	
	apply_tree_colors();
	tree_normal_vectors = [];
}

function generate_tree_root(tree_center, tree_height, tree_radius, tree_root_height){
		var shape_top = [];
		var shape_bot = [];
		var inc = 360 / num_tree_faces_root;
		
		var r = tree_radius / 3;
		
		for(let i = 0; i < num_tree_faces_root; ++i){
			shape_top.push( vec4( tree_center[0] + ( r * Math.cos(i * inc * deg_to_rad) ) + tree_bump_random(), tree_center[1] + tree_root_height, tree_center[2] + ( r * Math.sin(i * inc * deg_to_rad) ) + tree_bump_random(), 1.0 ) );
			shape_bot.push( vec4( tree_center[0] + ( r * Math.cos(i * inc * deg_to_rad) ) + tree_bump_random(), tree_center[1],					 tree_center[2] + ( r * Math.sin(i * inc * deg_to_rad) ) + tree_bump_random(), 1.0 ) );
		}
		
		for(let i = 0; i < shape_top.length; ++i){
			pointsArray.push(shape_top[i]);
		}
		
		let n = vec3(0.0, -1.0, 0.0);
		let l = vec3( - light[0], tree_center[1] - light[1], - light[2]);
				
		// find the magnitudes
		let n_mag = Math.sqrt( (n[0] * n[0]) + (n[1] * n[1]) + (n[2] * n[2]) );
		let l_mag = Math.sqrt( (l[0] * l[0]) + (l[1] * l[1]) + (l[2] * l[2]) );
		
		// find β
		let beta = ( (n[0] * l[0]) + (n[1] * l[1]) + (n[2] * l[2]) ) / (n_mag * l_mag);

		if(beta < 0){
			beta = 0;
		}
		
	for(let i = 0; i < shape_top.length; ++i){
		colorsArray.push( vec4(beta * tree_root_red, beta * tree_root_green, beta * tree_root_blue, 1.0) );
	}

	for(let i = 0; i < shape_bot.length; ++i){
		// get the points for the triangle
	
		let a = vec3(0.0, 0.0, 0.0);
		let b = shape_bot[(i + 1) % shape_bot.length];
		let c = shape_bot[i % shape_bot.length];
		
		// find the vectors
		let u = vec3((a[0] - b[0]), (a[1] - b[1]), (a[2] - b[2]));
		let v = vec3((c[0] - a[0]), (c[1] - a[1]), (c[2] - a[2]));
		
		let n = vec3( (u[1] * v[2]) - (u[2] * v[1]), (u[2] * v[0]) - (u[0] * v[2]),
		(u[0] * v[1]) - (u[1] * v[0]) );
		
		// now that we have the normal vector, find the light vector for this point
		let l = vec3( - light[0] + c[0], - light[1] + c[1], - light[2] + c[2]);
			
		// find the magnitudes
		let n_mag = Math.sqrt( (n[0] * n[0]) + (n[1] * n[1]) + (n[2] * n[2]) );
		let l_mag = Math.sqrt( (l[0] * l[0]) + (l[1] * l[1]) + (l[2] * l[2]) );
		
		// find β
		let beta = ( (n[0] * l[0]) + (n[1] * l[1]) + (n[2] * l[2]) ) / (n_mag * l_mag);
	
		if(beta < 0){
			beta = 0;
		}
		
		colorsArray.push( vec4(beta * tree_root_red, beta * tree_root_green, beta * tree_root_blue, 1.0) );
	}

		
		for(let i = 0; i < shape_bot.length; ++i){
			pointsArray.push(shape_bot[i]);
		}
						
		for(let i = 0; i <= shape_top.length; ++i){
			pointsArray.push( shape_top[i % shape_top.length] );
			pointsArray.push( shape_bot[i % shape_top.length] );
			pointsArray.push( shape_top[(i + 1) % shape_top.length] );
			pointsArray.push( shape_bot[(i + 1) % shape_bot.length] );

		// find the normal vector of the plane
		var a = shape_top[i % shape_top.length];
		var b = shape_bot[i % shape_top.length];
		var c = shape_bot[(i + 1) % shape_top.length];
		var d = shape_top[(i + 1) % shape_top.length]
		
		let u = vec3( (b[0] - a[0]), (b[1] - a[1]), (b[2] - a[2]) ); // b - a;
		let v = vec3( (c[0] - a[0]), (c[1] - a[1]), (c[2] - a[2]) ); // c - a;
		
		// find the normal vector
		let n = vec3( (u[1] * v[2]) - (u[2] * v[1]), (u[2] * v[0]) - (u[0] * v[2]),
		(u[0] * v[1]) - (u[1] * v[0]) );
				
		let l = light;
		// now that we have the normal vector, find the light vector for this point
		let lA = vec3( (a[0] - l[0]) , (a[1] - l[1]), (a[2] - l[2]) );		
		let lB = vec3( (b[0] - l[0]) , (b[1] - l[1]), (b[2] - l[2]) );		
		let lC = vec3( (c[0] - l[0]) , (c[1] - l[1]), (c[2] - l[2]) );		
		let lD = vec3( (d[0] - l[0]) , (d[1] - l[1]), (d[2] - l[2]) );		
		
		// find the magnitudes
		let n_mag = Math.sqrt( (n[0] * n[0]) + (n[1] * n[1]) + (n[2] * n[2]) );
		let lA_mag = Math.sqrt( (lA[0] * lA[0]) + (lA[1] * lA[1]) + (lA[2] * lA[2]) );
		let lB_mag = Math.sqrt( (lB[0] * lB[0]) + (lB[1] * lB[1]) + (lB[2] * lB[2]) );
		let lC_mag = Math.sqrt( (lC[0] * lC[0]) + (lC[1] * lC[1]) + (lC[2] * lC[2]) );
		let lD_mag = Math.sqrt( (lD[0] * lD[0]) + (lD[1] * lD[1]) + (lD[2] * lD[2]) );

		// find β
		let betaA = -( (n[0] * lA[0]) + (n[1] * lA[1]) + (n[2] * lA[2]) ) / (n_mag * lA_mag);
		let betaB = -( (n[0] * lB[0]) + (n[1] * lB[1]) + (n[2] * lB[2]) ) / (n_mag * lB_mag);
		let betaC = -( (n[0] * lC[0]) + (n[1] * lC[1]) + (n[2] * lC[2]) ) / (n_mag * lC_mag);
		let betaD = -( (n[0] * lD[0]) + (n[1] * lD[1]) + (n[2] * lD[2]) ) / (n_mag * lD_mag);
		
		betaA = Math.max(0, betaA);
		betaB = Math.max(0, betaB);
		betaC = Math.max(0, betaC);
		betaD = Math.max(0, betaD);
		
		colorsArray.push( vec4(betaA * tree_root_red, betaA * tree_root_green, betaA * tree_root_blue, 1.0) );
		colorsArray.push( vec4(betaB * tree_root_red, betaB * tree_root_green, betaB * tree_root_blue, 1.0) );
		colorsArray.push( vec4(betaC * tree_root_red, betaC * tree_root_green, betaC * tree_root_blue, 1.0) );
		colorsArray.push( vec4(betaD * tree_root_red, betaD * tree_root_green, betaD * tree_root_blue, 1.0) );
	}
}

function generate_unit_triangle_tree(a, b, c, d){
	var p1 = vec4(	a[0] / d, a[1] / d, a[2] / d, 1.0)
	var p2 = vec4(	b[0] / d, b[1] / d, b[2] / d, 1.0)
			
	let q = b[0] - p2[0];
	let w = b[1] - p2[1];
	let r = b[2] - p2[2];
	
	var diff = vec4(q, w, r, 1.0);
		
	var p3 = vec4(	c[0] / d, c[1] / d, c[2] / d, 1.0)
	var color_code = vec4(Math.random(), Math.random(), Math.random(), 1.0);

	generate_triangle_tree(p1, p2, p3, d, diff);
}

function generate_triangle_tree(a, b, c, d, diff){

	var cb_x_inc = c[0] - b[0];
	var cb_y_inc = c[1] - b[1];
	var cb_z_inc = c[2] - b[2];
	
	var ab_x_inc = a[0] - b[0];
	var ab_y_inc = a[1] - b[1];
	var ab_z_inc = a[2] - b[2];

	// start with the bottom
	let e = d;
	
	let q = diff[0];
	let w = diff[1];
	let r = diff[2];
	
	a = vec4( a[0] + q, a[1] + w, a[2] + r, 1.0 );
	b = vec4( b[0] + q, b[1] + w, b[2] + r, 1.0 );
	c = vec4( c[0] + q, c[1] + w, c[2] + r, 1.0 );
	
	var orig_a = a; // needed for normal vector calculation
	var orig_c = c; // needed for normal vector calculation
	
	var normal_a = a;
	var normal_b = b;
	var normal_c = c;

	var u;
	var v;

	for(let j = 0; j < d; ++j){
		pointsArray.push( vec4( tree_bump_random() + b[0],  tree_bump_random() + b[1],  tree_bump_random() + b[2], 1.0 ) );
		
		let cx = 0;
		let cz = 0;
		
		u = vec4( b[0] - a[0], b[1] - a[1], b[2] - a[2], 1.0 );
		v = vec4( c[0] - b[0], c[1] - b[1], c[2] - b[2], 1.0 );
	
		calc_normal_vector_tree(u, v); // three because we need all the points from the first triangle 
		calc_normal_vector_tree(u, v);
		calc_normal_vector_tree(u, v);
	
		var old_a;
		var old_c;
		
		for(let i = 0; i < e; ++i){
			normal_a = vec4(a[0] + cx + tree_bump_random(), a[1], a[2] + cz + tree_bump_random(), 1.0);
			normal_c = vec4(c[0] + cx + tree_bump_random(), c[1], c[2] + cz + tree_bump_random(), 1.0);
				
			pointsArray.push( normal_a );
			pointsArray.push( normal_c );
		
			if(first_triangle_in_tree_row_flag < 0){
				u = vec4( old_a[0] - normal_a[0], old_a[1] - normal_a[1], old_a[2] - normal_a[2], 1.0 );
				v = vec4( old_c[0] - old_a[0], old_c[1] - old_a[1], old_c[2] - old_a[2], 1.0 );
				
				calc_normal_vector_tree(u, v);
				
				u = vec4( old_c[0] - normal_a[0], old_c[1] - normal_a[1], old_c[2] - normal_a[2], 1.0 );
				v = vec4( normal_c[0] - old_c[0], normal_c[1] - old_c[1], normal_c[2] - old_c[2], 1.0 );
				
				calc_normal_vector_tree(u, v);				
				
			} else{
				first_triangle_in_tree_row_flag *= (-1);
			}
						
			cx += cb_x_inc;
			cz += cb_z_inc;
			
			old_a = normal_a;
			old_c = normal_c;
		}

		b = vec4( a[0], a[1], a[2], 1.0);
		c = vec4( b[0] + cb_x_inc, b[1] + cb_y_inc, b[2] + cb_z_inc, 1.0);		
		a = vec4( b[0] + ab_x_inc, b[1] + ab_y_inc, b[2] + ab_z_inc, 1.0);

		--e;
		first_triangle_in_tree_row_flag *= (-1);
	}
}

function apply_tree_colors(){
		for(let i = 0; i < tree_normal_vectors.length; ++i){
			var n = tree_normal_vectors[i];
			var p = pointsArray[i];
			
			var l = vec4( n[0] - p[0], n[1] - p[1], n[2] - p[2], 1.0 );
			
			var n_mag = Math.sqrt( (n[0] ** 2) + (n[1] ** 2) + (n[2] ** 2) );
			var l_mag = Math.sqrt( (l[0] ** 2) + (l[1] ** 2) + (l[2] ** 2) );
			
			var beta = ( ( n[0] * l[0] ) + ( n[1] * l[1] ) + ( n[2] * l[2] ) ) / (n_mag * l_mag);
			
			if(beta < 0){
				beta = 0;
			}
			
			colorsArray.push( vec4( tree_red * beta, tree_green * beta, tree_blue * beta, 1.0) );
		}
}

function tree_bump_random(){
	let min = -5;
	let max =  5;
	return (Math.random() * (max - min + 1) + min) / 1000;
}

function calc_normal_vector_tree(u, v){
	tree_normal_vectors.push( vec3( (u[1] * v[2]) - (u[2] * v[1]), (u[2] * v[0]) - (u[0] * v[2]),
		(u[0] * v[1]) - (u[1] * v[0]) ) );
}

var ship_center = vec3(0.8, 0.055, 0.8);
var ship_scale = 0.01;
var ship_top = 0.5;
var ship_bot = -0.75;

var pswr   = 134 / 255;
var pswg =  86 / 255;
var pswb  =  67 / 255;

var mast_height = 4 * ship_scale;
var mast_center_top = vec4( ship_center[0] + (0.5 * ship_scale), mast_height + ship_center[1] , ship_center[2], 1.0);
var mast_center_bot = vec4( ship_center[0] + (0.5 * ship_scale), ship_center[1] + ship_bot * ship_scale    , ship_center[2], 1.0);
var num_mast_sides = 20;

var sailr = 208 / 255;
var sailg = 198 / 255;
var sailb = 177 / 255;

var ship_normal_vectors = [];

var wood_points;
var ship_start;
var ship_rotation = 0;

var ship_shape_top = [
	vec4( ship_center[0] - (  ship_scale * 4 ),       ship_center[1] + ( ship_scale * ship_top ), ship_center[2] + ( ship_scale * 0 ),    1.0 ),
	vec4( ship_center[0] - (  ship_scale * (3.3) ),   ship_center[1] + ( ship_scale * ship_top ), ship_center[2] + ( ship_scale * -1.3 ),  1.0 ),
	vec4( ship_center[0] - (  ship_scale * (2) ),     ship_center[1] + ( ship_scale * ship_top ), ship_center[2] + ( ship_scale * -1.75 ),    1.0 ),
	vec4( ship_center[0] - (  ship_scale * (-1.5) ),  ship_center[1] + ( ship_scale * ship_top ), ship_center[2] + ( ship_scale * -2 ), 1.0 ),
	vec4( ship_center[0] - (  ship_scale * (-2) ),    ship_center[1] + ( ship_scale * ship_top ), ship_center[2] + ( ship_scale * 0 ),   1.0 ),
	vec4( ship_center[0] - (  ship_scale * (-1.5) ),  ship_center[1] + ( ship_scale * ship_top ), ship_center[2] + ( ship_scale * 2 ), 1.0 ),
	vec4( ship_center[0] - (  ship_scale * 2 ),       ship_center[1] + ( ship_scale * ship_top ), ship_center[2] + ( ship_scale * 1.75 ),    1.0 ),
	vec4( ship_center[0] - (  ship_scale * 3.3 ),     ship_center[1] + ( ship_scale * ship_top ), ship_center[2] + ( ship_scale * 1.3 ),  1.0 ),
	vec4( ship_center[0] - (  ship_scale * 4 ),       ship_center[1] + ( ship_scale * ship_top ), ship_center[2] + ( ship_scale * 0 ), 1.0 )
];

var ship_shape_bot = [
	vec4( ship_center[0] - (  ship_scale * 3 ),      ship_center[1] + ( ship_scale * ship_bot ), ship_center[2] + ( ship_scale * 0 ),      1.0 ),
	vec4( ship_center[0] - (  ship_scale * 2.475 ),  ship_center[1] + ( ship_scale * ship_bot ), ship_center[2] + ( ship_scale * -0.975 ),  1.0 ),
	vec4( ship_center[0] - (  ship_scale * 1.5 ),    ship_center[1] + ( ship_scale * ship_bot ), ship_center[2] + ( ship_scale * -1.313 ),    1.0 ),
	vec4( ship_center[0] - (  ship_scale * -1.125 ), ship_center[1] + ( ship_scale * ship_bot ), ship_center[2] + ( ship_scale * -1.5 ), 1.0 ),
	vec4( ship_center[0] - (  ship_scale * -1.5 ),   ship_center[1] + ( ship_scale * ship_bot ), ship_center[2] + ( ship_scale * 0 ),   1.0 ),
	vec4( ship_center[0] - (  ship_scale * -1.125 ), ship_center[1] + ( ship_scale * ship_bot ), ship_center[2] + ( ship_scale * 1.5 ), 1.0 ),
	vec4( ship_center[0] - (  ship_scale * 1.5 ),    ship_center[1] + ( ship_scale * ship_bot ), ship_center[2] + ( ship_scale * 1.313 ),    1.0 ),
	vec4( ship_center[0] - (  ship_scale * 2.475 ),  ship_center[1] + ( ship_scale * ship_bot ), ship_center[2] + ( ship_scale * 0.975 ),  1.0 ),
	vec4( ship_center[0] - (  ship_scale * 3 ),      ship_center[1] + ( ship_scale * ship_bot ), ship_center[2] + ( ship_scale * 0 ),      1.0 )
];


function generate_pirate_ship(){
	wood_points = 118;

	create_main_body();
	create_mast();

	sail_start = pointsArray.length;
	create_sail();
	
	apply_ship_colors_wood();
	apply_ship_colors_sails();
}

function create_sail(){
	var p1 = vec4( ship_center[0] -  ( 1 * ship_scale),  ship_center[1] + (1.5  * ship_scale) , ship_center[2] + (2 * ship_scale), 1.0 );
	var p2 = vec4( ship_center[0] -  ( 1 * ship_scale),  ship_center[1] + (4 * ship_scale),     ship_center[2] + (2 * ship_scale), 1.0 );
	var p3 = vec4( ship_center[0] +-  ( 1 * ship_scale),  ship_center[1] + (4 * ship_scale),     ship_center[2] + (-2 * ship_scale), 1.0 );
	var p4 = vec4( ship_center[0] + - ( 1 * ship_scale),  ship_center[1] + (1.5 * ship_scale),   ship_center[2] + (-2 * ship_scale), 1.0 );

	var p5 = vec4( ship_center[0] +  - (1.5 * ship_scale), ship_center[1] + 2.25 * ship_scale, ship_center[2] + (1.25 * ship_scale), 1.0 );
	var p6 = vec4( ship_center[0] +  - (1.5 * ship_scale), ship_center[1] + 3.25 * ship_scale, ship_center[2] + (1.25 * ship_scale), 1.0 );
	var p7 = vec4( ship_center[0] +  - (1.5 * ship_scale), ship_center[1] + 3.25 * ship_scale, ship_center[2] + (-1.25 * ship_scale), 1.0 );
	var p8 = vec4( ship_center[0] +  - (1.5 * ship_scale), ship_center[1] + 2.25 * ship_scale, ship_center[2] + (-1.25 * ship_scale), 1.0 );

/*
	var p1 = vec4( ship_center[0] +  ( 2 * ship_scale),  ship_center[1] + (1.5  * ship_scale) , ship_center[2] + 1 * ship_scale, 1.0 );
	var p2 = vec4( ship_center[0] +  ( 2 * ship_scale),  ship_center[1] + (4 * ship_scale),     ship_center[2] + 1 * ship_scale, 1.0 );
	var p3 = vec4( ship_center[0] +  (-2 * ship_scale),  ship_center[1] + (4 * ship_scale),     ship_center[2] + 1 * ship_scale, 1.0 );
	var p4 = vec4( ship_center[0] +  (-2 * ship_scale),  ship_center[1] + (1.5 * ship_scale),   ship_center[2] + 1 * ship_scale, 1.0 );

	var p5 = vec4( ship_center[0] +   1.25 * ship_scale, ship_center[1] + 2.25 * ship_scale, ship_center[2] + 1.5 * ship_scale, 1.0 );
	var p6 = vec4( ship_center[0] +   1.25 * ship_scale, ship_center[1] + 3.25 * ship_scale, ship_center[2] + 1.5 * ship_scale, 1.0 );
	var p7 = vec4( ship_center[0] +  -1.25 * ship_scale, ship_center[1] + 3.25 * ship_scale, ship_center[2] + 1.5 * ship_scale, 1.0 );
	var p8 = vec4( ship_center[0] +  -1.25 * ship_scale, ship_center[1] + 2.25 * ship_scale, ship_center[2] + 1.5 * ship_scale, 1.0 );
*/
	
	// push the points into the thing
	create_sail_polygon(p5, p6, p7, p8); // front

	create_sail_polygon(p2, p6, p5, p1);
	create_sail_polygon(p3, p7, p6, p2);
	create_sail_polygon(p4, p8, p7, p3);
	create_sail_polygon(p1, p5, p8, p4);
}

function create_sail_polygon(a, b, c, d){
			let u = vec4( (a[0] - c[0]), (a[1] - c[1]), (a[2] - c[2]), 1.0 );
			let v = vec4( (b[0] - a[0]), (b[1] - a[1]), (b[2] - a[2]), 1.0 );
			
			let n = get_normal_vector(u, v);

			ship_normal_vectors.push(n);
			ship_normal_vectors.push(n);
			ship_normal_vectors.push(n);
			ship_normal_vectors.push(n);
			
			pointsArray.push( a );
			pointsArray.push( b );
			pointsArray.push( c );
			pointsArray.push( d );
}

function create_mast(){
	// mast_base = ship_top;
	var mast_radius = 0.3 * ship_scale;
	mast_height = 6 * ship_scale;
	
	// create a cylinder of height ( mast_height - ship_top ) and radius ( 0.3 * scale)
	var s_top = [];
	var s_bot = [];
	var mast_inc = 360 / num_mast_sides;
	
	for(let i = 0; i < num_mast_sides; ++i){
		let t = vec4(
		mast_center_top[0] + mast_radius * Math.cos(i * mast_inc * deg_to_rad),
		mast_center_top[1],
		mast_center_top[2] + mast_radius * Math.sin(i * mast_inc * deg_to_rad),
		1.0		);
		
		let q = vec4( 
		mast_center_bot[0] + mast_radius * Math.cos(i * mast_inc * deg_to_rad),
		mast_center_bot[1],
		mast_center_bot[2] + mast_radius * Math.sin(i * mast_inc * deg_to_rad),
		1.0 );
				
		s_top.push( t );
		s_bot.push( q );

	}

	// push the shapes and get the normal_vectors
			// push the top
		for(let i = 0; i < s_top.length; ++i){
			pointsArray.push( s_top[i] );	
			ship_normal_vectors.push( vec3(0.0, 1.0, 0.0) );
		}
		
		// push the bottom
		for(let i = 0; i < s_bot.length; ++i){
			pointsArray.push( s_bot[i] );
			ship_normal_vectors.push( vec3(0.0, -1.0, 0.0) );
		}
		
		// push the side
		for(let i = 0; i <= s_top.length; ++i){
			let a = s_top[i % s_top.length];
			let b = s_bot[i % s_top.length];
			let c = s_top[ (i + 1) % s_top.length ];
			
			let u = vec4( (b[0] - a[0]), (b[1] - a[1]), (b[2] - a[2]), 1.0 );
			let v = vec4( (a[0] - c[0]), (a[1] - c[1]), (a[2] - c[2]), 1.0 );
			
			let n = get_normal_vector(v, u);

			ship_normal_vectors.push(n);
			ship_normal_vectors.push(n);
			
			pointsArray.push( a );
			pointsArray.push( b );
		}
}

function create_main_body(){
		// push the top
		for(let i = 0; i < ship_shape_top.length; ++i){
			pointsArray.push( ship_shape_top[i] );	
			ship_normal_vectors.push( vec3(0.0, 1.0, 0.0) );
		}
		
		// push the bottom
		for(let i = 0; i < ship_shape_bot.length; ++i){
			pointsArray.push( ship_shape_bot[i] );
			ship_normal_vectors.push( vec3(0.0, -1.0, 0.0) );
		}
		
		// push the side
		for(let i = 0; i < ship_shape_top.length; ++i){
			let a = ship_shape_top[i];
			let b = ship_shape_bot[i];
			let c = ship_shape_top[ (i + 1) % ship_shape_top.length ];
			
			let u = vec4( (b[0] - a[0]), (b[1] - a[1]), (b[2] - a[2]), 1.0 );
			let v = vec4( (c[0] - a[0]), (c[1] - a[1]), (c[2] - a[2]), 1.0 );
			
			let n = get_normal_vector(u, v);

			ship_normal_vectors.push(n);
			ship_normal_vectors.push(n);
			
			pointsArray.push( ship_shape_top[i] );
			pointsArray.push( ship_shape_bot[i] );
		}
}

var sail_start;

function apply_ship_colors_sails(){
	
	for(let i = 0; i < 20; ++i){
		let p = pointsArray[sail_start + i];
		let n = ship_normal_vectors[i];

		// get the vector from the point to the light
		let l = vec4( light[0] - p[0], light[1] - p[1], light[2] - p[2], 1.0);

		// find the magnitudes
		let nm = Math.sqrt( ( n[0] * n[0] ) + ( n[1] * n[1] ) + ( n[2] * n[2] ) );
		let lm = Math.sqrt( ( l[0] * l[0] ) + ( l[1] * l[1] ) + ( l[2] * l[2] ) );

		// find beta
		
		let beta = ( ( n[0] * l[0] ) + ( n[1] * l[1] ) + ( n[2] * l[2] ) ) / ( nm * lm );
		if(beta < 0){
			beta *= (-1);
		}
		beta = Math.max(beta, 0);
		
		colorsArray.push( vec4( beta * sailr, beta * sailg, beta * sailb, 1.0 ) );
	}
}

function apply_ship_colors_wood(){
	//wood_points
	for(let i = 0; i < 118; ++i){
		let p = pointsArray[i + ship_start];
		let n = ship_normal_vectors[i];

		// get the vector from the point to the light
		let l = vec4( light[0] - p[0], light[1] - p[1], light[2] - p[2], 1.0);
		
		// find the magnitudes
		let nm = Math.sqrt( ( n[0] * n[0] ) + ( n[1] * n[1] ) + ( n[2] * n[2] ) );
		let lm = Math.sqrt( ( l[0] * l[0] ) + ( l[1] * l[1] ) + ( l[2] * l[2] ) );

		// find beta
		let beta = ( ( n[0] * l[0] ) + ( n[1] * l[1] ) + ( n[2] * l[2] ) ) / ( nm * lm );
		beta = Math.max(beta, 0);
		
		colorsArray.push( vec4( beta * pswr, beta * pswg, beta * pswb, 1.0 ) );
	}
}


function get_normal_vector(u, v){
	return vec4( 
	- ((u[1] * v[2]) - (u[2] * v[1])),
	- ((u[2] * v[0]) - (u[0] * v[2])),
	- ((u[0] * v[1]) - (u[1] * v[0])),
	1.0
	);	
}

var ship_theta = 0;
var ship_end;


function update_ship(){
	var ship = [];
	
	for(let i = ship_start; i < ship_end; ++i){
		let p = pointsArray[i];
		let point = vec4(p[0], p[1] + Math.cos( water_trig[ ship_theta ] ) / 150, p[2], 1.0);
		
		
		ship.push(point);
	}
	
	ship_theta = ( ship_theta + 1 ) % water_trig.length;
	
	gl.bindBuffer(gl.ARRAY_BUFFER, vBuffer);
    gl.bufferSubData(gl.ARRAY_BUFFER, (ship_start * 16), flatten(ship));
}








var vBuffer;

window.onload = function init() {
    canvas = document.getElementById("gl-canvas");

    gl = WebGLUtils.setupWebGL(canvas);
    if (!gl) {
        alert("WebGL isn't available");
    }

    gl.viewport(0, 0, canvas.width, canvas.height);

    aspect = canvas.width / canvas.height;

    gl.clearColor(1.0, 1.0, 1.0, 1.0);

    gl.enable(gl.DEPTH_TEST);

    //
    //  Load shaders and initialize attribute buffers
    //
    var program = initShaders(gl, "vertex-shader", "fragment-shader");
    gl.useProgram(program);

	generate_trig();

    genTerrainData(terrRows, terrCols);
    prepMesh(terrRows, terrCols);

    colorCube();
	
	grass_start = pointsArray.length;
	generate_grass_underground();
		
    mountain_start = pointsArray.length;
    generate_mountain();

	water_start = pointsArray.length;
	generate_water_first_time();
	water_length = pointsArray.length - water_start; 

	sand_start = pointsArray.length;
	generate_sand();	

	tree1_start = pointsArray.length;
	generate_tree(tree_center1, tree_height1, tree_radius1, tree_root_height1);	
	
	tree2_start = pointsArray.length;
	generate_tree(tree_center2, tree_height2, tree_radius2, tree_root_height2);

	tree3_start = pointsArray.length;
	generate_tree(tree_center3, tree_height3, tree_radius3, tree_root_height3);
	
	tree4_start = pointsArray.length;
	generate_tree(tree_center4, tree_height4, tree_radius4, tree_root_height4);

	tree5_start = pointsArray.length;
	generate_tree(tree_center5, tree_height5, tree_radius5, tree_root_height5);

	tree6_start = pointsArray.length;
	generate_tree(tree_center6, tree_height6, tree_radius6, tree_root_height6);

	ship_start = pointsArray.length; 
    generate_pirate_ship();

	ship_end = pointsArray.length;

    var cBuffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, cBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, flatten(colorsArray), gl.STATIC_DRAW);

    var vColor = gl.getAttribLocation(program, "vColor");
    gl.vertexAttribPointer(vColor, 4, gl.FLOAT, false, 0, 0);
    gl.enableVertexAttribArray(vColor);

    vBuffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, vBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, flatten(pointsArray), gl.STATIC_DRAW);

    var vPosition = gl.getAttribLocation(program, "vPosition");
    gl.vertexAttribPointer(vPosition, 4, gl.FLOAT, false, 0, 0);
    gl.enableVertexAttribArray(vPosition);

    modelView = gl.getUniformLocation(program, "modelView");
    projection = gl.getUniformLocation(program, "projection");

    document.getElementById("Button9").onclick = function() {
        pause_flag *= (-1);
    };
//	document.getElementById("Button10").onclick = function() {
//        theta_inc *= (-1);
//    };
    render();
}

var theta_inc = Math.PI / 450;

var fovy = 45.0; // Field-of-view in Y direction angle (in degrees)
var aspect; // Viewport aspect ratio

var eye_theta_for_test = 0;

var mountain_cam_radius = 0.45;
var mountain_phi = pi;
var mountain_cam_height = 1 / 5;
var fluff = 1;

var t = 0.0006;

var mountain_cam_flag = 1;

// initialize the eyes and 
var bee_const = 0; // maybe I'll use, maybe I won't, depends on the bee

var mountain_eye_x =  mountain_center[0] + (mountain_cam_radius * Math.sin(mountain_phi));
var mountain_at_x = mountain_center[0] + ( mountain_cam_radius - t ) * Math.sin( ( mountain_phi + (2 * theta_inc) ));
var mountain_eye_y = mountain_cam_height;
var mountain_at_y = mountain_cam_height + (1 / 1000);
var mountain_eye_z = mountain_center[2] + mountain_cam_radius * Math.cos(mountain_phi);
var mountain_at_z = mountain_center[2] + ( mountain_cam_radius - t) * Math.cos(mountain_phi + (2 * theta_inc) );

// position: ( -0.5, 0.2, -0.95 )
// end: (-0.7265057831928067, 0.47500000000000026, -0.5399390808633996)


var mountain_animation_frame = 1;
var first_mountain_flag = 21;

var z_check_flag = -0.6701981033940475;

var island_lap_height = 0.7;

var island_fly_theta = (pi / 2) - (pi / 35);
var island_fly_radius = ( 7 * Math.sin( 2 * island_fly_theta ) / 3 ) - 1;

var speed = 1;

var final_x_stretch = 1 / 3;
var final_z_stretch = -0.47;

var render = function() {
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    // change the camera location or position

//	eye = vec3(mountain_eye_x, mountain_eye_y, mountain_eye_z );
//	at = vec3(mountain_at_x , mountain_at_y,  mountain_at_z );
	

    if (pause_flag > 0) {
		// fly around the mountain
		if( mountain_phi  + ( 30 * theta_inc ) <= ( 7 * pi / 2 ) ) {
			mountain_cam_height += speed * 1 / 4000;
			mountain_cam_radius -= speed * 0.0002;
			mountain_phi += speed * theta_inc;

			mountain_eye_x = mountain_center[0] + ( mountain_cam_radius * Math.sin(mountain_phi) );
			mountain_eye_y = mountain_cam_height;
			mountain_eye_z = mountain_center[2] + ( mountain_cam_radius * Math.cos(mountain_phi) );

			mountain_at_x  = mountain_center[0] + ( mountain_cam_radius - t ) * Math.sin( ( mountain_phi + (2 * theta_inc) ) );
			mountain_at_y  = mountain_cam_height;
			mountain_at_z  = mountain_center[2] + ( mountain_cam_radius - t ) * Math.cos( ( mountain_phi + (2 * theta_inc) ) );

			eye = vec3(mountain_eye_x, mountain_eye_y, mountain_eye_z );
			at = vec3(mountain_at_x , mountain_at_y,  mountain_at_z );
		} else{
				// do the lap around the isladnd
				if(island_lap_height >= 1 / 5){
					eye = vec3(( island_fly_radius * Math.sin(island_fly_theta) ), island_lap_height, island_fly_radius * Math.cos( island_fly_theta ) );
					at = vec3(( island_fly_radius * Math.sin(island_fly_theta + (2 * theta_inc) ) ), island_lap_height - (4 / 1000), island_fly_radius * Math.cos( island_fly_theta + (2 * theta_inc) ) );
					
					island_lap_height -= 1 / ( 1200 );
					island_fly_theta += theta_inc;
					phi -= theta_inc;

				} else{
					if(final_x_stretch > -0.4 || final_z_stretch > -0.8){
						eye = vec3(final_x_stretch, 0.2, final_z_stretch);
						at = vec3(-1.0, 0.2, -1.0);							
						
						let o = 10;

						final_x_stretch -= (0.04 / o);
						final_z_stretch -= (0.025 / o);						
					} else{
						// reset the terrain flags
						mountain_cam_height = 1 / 5;
						mountain_cam_radius = 0.45;
						mountain_phi = pi;
						z_check_flag = -0.6701981033940475;
						island_lap_height = 0.7;
						phi = 270.0 * Math.PI / 180.0;						
						island_fly_theta = ( pi / 2 ) - ( pi / 35 );
						island_fly_radius = ( 7 * Math.sin( 2 * island_fly_theta ) / 3 ) - 1;
						final_x_stretch = 1 / 3;
						final_z_stretch = -0.47;
					}
				}
		}
    }

	mvMatrix = lookAt(eye, at, up);
    pMatrix = perspective(fovy, aspect, far, near); // symmetrical around z

    /*recall P' = P*MV*p*/

    // pass info to the shader
    gl.uniformMatrix4fv(modelView, false, flatten(mvMatrix));
    gl.uniformMatrix4fv(projection, false, flatten(pMatrix));
  
    for (var i = 0; i < terrRows; ++i) {
    	gl.drawArrays(gl.LINE_STRIP, i*terrCols, terrCols);
    }
    for (var i = 0; i < terrCols; ++i) {
    	gl.drawArrays(gl.LINE_STRIP, i*terrRows+(index/2), terrRows);
    }
   
    gl.drawArrays(gl.TRIANGLES, sky_start, numVertices);

	gl.drawArrays(gl.TRIANGLE_STRIP, grass_start, 4);
	
    var point_counter = mountain_start;

    for (let i = 0; i < sides; ++i) {
        for (let j = num_rows_per_triangle; j > 0; --j) {
            gl.drawArrays(gl.TRIANGLE_STRIP, point_counter, (2 * j) + 1);
            point_counter += (2 * j) + 1;
        }
    }

	let water_count = water_start;
	
	for(let i = 0; i < rows_of_water; ++i){
		gl.drawArrays(gl.TRIANGLE_STRIP, water_count, water_length / rows_of_water);		
		water_count +=  water_length / rows_of_water;
	}

	let sand_count = sand_start;
	for(let i = 0; i < rows_of_sand; ++i){
		gl.drawArrays(gl.TRIANGLE_STRIP, sand_count, sand_length / rows_of_sand);		
		sand_count += sand_length / rows_of_sand;
	}
	
	let tree_inc = (2 * num_tree_faces_root) + ( 4 * ( num_tree_faces_root + 1 ) );

	gl.drawArrays( gl.TRIANGLE_FAN, tree1_start, num_tree_faces_root );
	gl.drawArrays( gl.TRIANGLE_FAN, tree1_start + num_tree_faces_root, num_tree_faces_root );
	gl.drawArrays( gl.TRIANGLE_STRIP, tree1_start + ( 2 * num_tree_faces_root ), 4 * num_tree_faces_root);

	for(let i = 0; i < sides_per_tree; ++i){
		for(let j = num_rows_per_tree_triangle; j > 0; --j){
			gl.drawArrays(gl.TRIANGLE_STRIP, tree1_start + tree_inc, (2 * j) + 1);
			tree_inc += (2 * j) + 1;
		}
	}

	tree_inc = (2 * num_tree_faces_root) + ( 4 * ( num_tree_faces_root + 1 ) );

	gl.drawArrays( gl.TRIANGLE_FAN, tree2_start, num_tree_faces_root );
	gl.drawArrays( gl.TRIANGLE_FAN, tree2_start + num_tree_faces_root, num_tree_faces_root );
	gl.drawArrays( gl.TRIANGLE_STRIP, tree2_start + ( 2 * num_tree_faces_root ), 4 * num_tree_faces_root);

	for(let i = 0; i < sides_per_tree; ++i){
		for(let j = num_rows_per_tree_triangle; j > 0; --j){
			gl.drawArrays(gl.TRIANGLE_STRIP, tree2_start + tree_inc, (2 * j) + 1);
			tree_inc += (2 * j) + 1;
		}
	}

	tree_inc = (2 * num_tree_faces_root) + ( 4 * ( num_tree_faces_root + 1 ) );

	gl.drawArrays( gl.TRIANGLE_FAN, tree3_start, num_tree_faces_root );
	gl.drawArrays( gl.TRIANGLE_FAN, tree3_start + num_tree_faces_root, num_tree_faces_root );
	gl.drawArrays( gl.TRIANGLE_STRIP, tree3_start + ( 2 * num_tree_faces_root ), 4 * num_tree_faces_root);

	for(let i = 0; i < sides_per_tree; ++i){
		for(let j = num_rows_per_tree_triangle; j > 0; --j){
			gl.drawArrays(gl.TRIANGLE_STRIP, tree3_start + tree_inc, (2 * j) + 1);
			tree_inc += (2 * j) + 1;
		}
	}
	
		tree_inc = (2 * num_tree_faces_root) + ( 4 * ( num_tree_faces_root + 1 ) );

	gl.drawArrays( gl.TRIANGLE_FAN, tree4_start, num_tree_faces_root );
	gl.drawArrays( gl.TRIANGLE_FAN, tree4_start + num_tree_faces_root, num_tree_faces_root );
	gl.drawArrays( gl.TRIANGLE_STRIP, tree4_start + ( 2 * num_tree_faces_root ), 4 * num_tree_faces_root);

	for(let i = 0; i < sides_per_tree; ++i){
		for(let j = num_rows_per_tree_triangle; j > 0; --j){
			gl.drawArrays(gl.TRIANGLE_STRIP, tree4_start + tree_inc, (2 * j) + 1);
			tree_inc += (2 * j) + 1;
		}
	}

	tree_inc = (2 * num_tree_faces_root) + ( 4 * ( num_tree_faces_root + 1 ) );

	gl.drawArrays( gl.TRIANGLE_FAN, tree5_start, num_tree_faces_root );
	gl.drawArrays( gl.TRIANGLE_FAN, tree5_start + num_tree_faces_root, num_tree_faces_root );
	gl.drawArrays( gl.TRIANGLE_STRIP, tree5_start + ( 2 * num_tree_faces_root ), 4 * num_tree_faces_root);

	for(let i = 0; i < sides_per_tree; ++i){
		for(let j = num_rows_per_tree_triangle; j > 0; --j){
			gl.drawArrays(gl.TRIANGLE_STRIP, tree5_start + tree_inc, (2 * j) + 1);
			tree_inc += (2 * j) + 1;
		}
	}

	tree_inc = (2 * num_tree_faces_root) + ( 4 * ( num_tree_faces_root + 1 ) );

	gl.drawArrays( gl.TRIANGLE_FAN, tree6_start, num_tree_faces_root );
	gl.drawArrays( gl.TRIANGLE_FAN, tree6_start + num_tree_faces_root, num_tree_faces_root );
	gl.drawArrays( gl.TRIANGLE_STRIP, tree6_start + ( 2 * num_tree_faces_root ), 4 * num_tree_faces_root);

	for(let i = 0; i < sides_per_tree; ++i){
		for(let j = num_rows_per_tree_triangle; j > 0; --j){
			gl.drawArrays(gl.TRIANGLE_STRIP, tree6_start + tree_inc, (2 * j) + 1);
			tree_inc += (2 * j) + 1;
		}
	}

    gl.drawArrays( gl.TRIANGLE_FAN, ship_start, 9 ); // ship_top
    gl.drawArrays( gl.TRIANGLE_FAN, ship_start + 9, 9 ); // ship_bot
	gl.drawArrays( gl.TRIANGLE_STRIP, ship_start + 18, 18 ); // ship_sides
		
	let mast_start = ship_start + 36;
	for(let i = 0; i < num_mast_sides; ++i){
		gl.drawArrays( gl.TRIANGLE_STRIP, mast_start , 6);
		mast_start += 4;
    }

	let sail_count = sail_start;
	for(let i = 0; i < 5; ++i){
		gl.drawArrays( gl.TRIANGLE_FAN, sail_count, 4);
		sail_count += 4;
	}

	generate_water();
	update_ship();
	
    requestAnimFrame(render);
}
