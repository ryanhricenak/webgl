"use strict";

var canvas;
var gl;

var points = [];
var colors = [];

var axis = 0;
var theta = [0, 0, 0];
var thetaLoc;

var r = 0.1; // this is the radius for the shape
var hr = Math.sqrt(2) * (r / 2);

var cube_theta = 0;
var axis_theta = 0;
var deg_to_rad = Math.PI / 180;

var orig_center = [];
var center = [];

var num_to_start = 20; // needed to first initialize the buffer
var num_of_cubes = 20; 

var old_num = 20;
var flag = 1;

var vBuffer;
var cBuffer;

var original_centers = [];
var centers = [];
var color_codes = [];

// from here on out, it's survival time

var near = -0.5;
var far = 0.5;
var eye_radius = 1.0;
var eye_theta  = 45.0 * Math.PI/180.0;
var phi    = 0.0 * Math.PI/180.0;
var dr = 5.0 * Math.PI/180.0;

var left = -1.0;
var right = 1.0;
var ytop = 1.0;
var bottom = -1.0;

var mvMatrix, pMatrix;
var modelView, projection;
var eye;

const at = vec3(0.0, 0.0, 0.0);/*look at origin*/
const up = vec3(0.0, 1.0, 0.0);/*y-axis direction as up*/

var  fovy = 45.0;  // Field-of-view in Y direction angle (in degrees)
var  aspect;       // Viewport aspect ratio


window.onload = function init()
{
    canvas = document.getElementById( "gl-canvas" );

    gl = WebGLUtils.setupWebGL( canvas );
    if ( !gl ) { alert( "WebGL isn't available" ); }
    
	aspect = canvas.width / canvas.height;
	
    gl.viewport( 0, 0, canvas.width, canvas.height );
    gl.clearColor( 0.0, 0.0, 0.0, 1.0 );

	for(let i = 0; i < num_of_cubes; ++i){
		let v = vec3(gr() , yr(), gr());
		original_centers.push(v);
		centers.push( v );			
	}
	
	create_cube();
	
    gl.enable(gl.DEPTH_TEST);

    //
    //  Load shaders and initialize attribute buffers
    //
    var program = initShaders( gl, "vertex-shader", "fragment-shader" );
    gl.useProgram( program );

    cBuffer = gl.createBuffer();
    gl.bindBuffer( gl.ARRAY_BUFFER, cBuffer );
    gl.bufferData( gl.ARRAY_BUFFER, flatten(colors), gl.STATIC_DRAW );

    var vColor = gl.getAttribLocation( program, "vColor" );
    gl.vertexAttribPointer( vColor, 4, gl.FLOAT, false, 0, 0 );
    gl.enableVertexAttribArray( vColor );

    thetaLoc = gl.getUniformLocation(program, "theta");


    vBuffer = gl.createBuffer();
    gl.bindBuffer( gl.ARRAY_BUFFER, vBuffer );
    gl.bufferData( gl.ARRAY_BUFFER, flatten(points), gl.STATIC_DRAW );


    var vPosition = gl.getAttribLocation( program, "vPosition" );
    gl.vertexAttribPointer( vPosition, 4, gl.FLOAT, false, 0, 0 );
    gl.enableVertexAttribArray( vPosition );



	modelView = gl.getUniformLocation( program, "modelView");
	projection = gl.getUniformLocation( program, "projection" );







    var m = document.getElementById("num_cubes");
		m.addEventListener("click", function(){
     	num_of_cubes = m.value;
		
		if(old_num === num_of_cubes || flag === 1){
			flag = 0;
		} else{
			old_num = num_of_cubes;
			add_centers();
		}		
    });

	num_of_cubes = 1;
    render();
}

// function that adds cubes to the canvas
function add_centers(){
	centers = [];
	original_centers = [];
	colors = [];
	
	for(let i = 0; i < num_of_cubes; ++i){
			let v = vec3(gr(), yr(), gr());
			centers.push(v);
			original_centers.push(v);
	}
	create_cube2();
}

// function that returns a random value between -0.65 and 0.65
// used for x and z values only
function gr(){
	let max = 0.65;
	let min = -0.65;
	return Math.random() * (max - min) + min;
}

// function that returns a random value between -0.9 and 0.9
// used for the y values only
function yr(){
	let max = 0.9;
	let min = -0.9;
	return Math.random() * (max - min) + min;
}

// function to first initialize the list
function create_cube(){
	
	points = [];
	
	for(let i = 0; i < 20 * 16; ++i){

		let red = Math.random();
		let gre = Math.random();
		let blu = Math.random();
		let color_code = vec4(red, gre, blu, 1.0);

		for(let j = 0; j < 16; ++j){
			points.push( vec4(0.0, 0.0, 0.0, 0.0));
			colors.push(color_code);			
		}
	}
}

function create_cube2(){

	var p = [];
	
	for(let i = 0; i < num_of_cubes; ++i){
		let a = centers[i][0];
		let b = centers[i][1];
		let c = centers[i][2];
	
		// points for the cube that is centered around center
		let p1 = vec4(a + r * Math.cos(  cube_theta * deg_to_rad) ,			 b + hr, c + r * Math.sin(  cube_theta * deg_to_rad), 1.0);
		let p2 = vec4(a + r * Math.cos( (cube_theta + 90) * deg_to_rad) ,	 b + hr, c + r * Math.sin( (cube_theta + 90) * deg_to_rad), 1.0);
		let p3 = vec4(a + r * Math.cos( (cube_theta + 180) * deg_to_rad) ,	 b + hr, c + r * Math.sin( (cube_theta + 180) * deg_to_rad), 1.0);
		let p4 = vec4(a + r * Math.cos( (cube_theta + 270) * deg_to_rad) ,	 b + hr, c + r * Math.sin( (cube_theta + 270) * deg_to_rad), 1.0);
		let p5 = vec4(a + r * Math.cos(cube_theta * deg_to_rad) ,			 b - hr, c + r * Math.sin(  cube_theta * deg_to_rad), 1.0);
		let p6 = vec4(a + r * Math.cos( (cube_theta + 90) * deg_to_rad) , 	 b - hr, c + r * Math.sin( (cube_theta + 90) * deg_to_rad), 1.0);
		let p7 = vec4(a + r * Math.cos( (cube_theta + 180) * deg_to_rad) , 	 b - hr, c + r * Math.sin( (cube_theta + 180) * deg_to_rad), 1.0);
		let p8 = vec4(a + r * Math.cos( (cube_theta + 270) * deg_to_rad) , 	 b - hr, c + r * Math.sin( (cube_theta + 270) * deg_to_rad), 1.0);	

		// push the needed points
		p.push(p1);
		p.push(p2);
		p.push(p3);
		p.push(p4);

		p.push(p5);
		p.push(p6);
		p.push(p7);
		p.push(p8);

		p.push(p1);
		p.push(p5);
	
		p.push(p2);
		p.push(p6);
	
		p.push(p3);
		p.push(p7);
	
		p.push(p4);
		p.push(p8);
	
		//for(let t = 0; t < 16; ++t){
			//colors.push( vec4(0.0, 0.0, 1.0, 1.0) );
		//}
	
	}

	gl.bindBuffer(gl.ARRAY_BUFFER, vBuffer);
	gl.bufferSubData(gl.ARRAY_BUFFER, 0, flatten(p));
}

function render()
{
    gl.clear( gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

	eye = vec3( 1.5, 1.5, 1.5 );
	mvMatrix = lookAt(eye, at, up);
	pMatrix = perspective(fovy, aspect, far, near);

	gl.uniformMatrix4fv (modelView, false, flatten(mvMatrix));
	gl.uniformMatrix4fv (projection, false, flatten(pMatrix));


	theta[1] += Math.PI / 5; 
    gl.uniform3fv(thetaLoc, theta);

	cube_theta += 2.5;
	
	let cube_offset = 0;
	for(let i = 0; i < num_of_cubes; ++i){
		gl.drawArrays( gl.LINE_LOOP, cube_offset + 0, 4 );
		gl.drawArrays( gl.LINE_LOOP, cube_offset + 4, 4 );
		gl.drawArrays( gl.LINES, cube_offset + 8, 8 );
		cube_offset += 16;

		create_cube2();
	}

    requestAnimFrame( render );
}
