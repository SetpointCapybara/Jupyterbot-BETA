import { Object3D, Vector3, BoxBufferGeometry, Color, Mesh, MeshBasicMaterial, PerspectiveCamera,
Scene, WebGLRenderer, AmbientLight, DirectionalLight, HemisphereLight, MeshStandardMaterial,
AxesHelper, GridHelper, Matrix4, SphereBufferGeometry, CylinderBufferGeometry, Group, LoadingManager, MeshPhysicalMaterial, Vector2, FrontSide,
	BackSide, DoubleSide, PMREMGenerator, TextureLoader, PointLight, UVMapping, CubeReflectionMapping, CubeRefractionMapping,
	EquirectangularReflectionMapping, EquirectangularRefractionMapping, CubeUVReflectionMapping,
	CubeUVRefractionMapping, RepeatWrapping, ClampToEdgeWrapping, MirroredRepeatWrapping,
	NearestFilter, LinearFilter, NearestMipmapNearestFilter, NearestMipmapLinearFilter, LinearMipmapNearestFilter,
	LinearMipmapLinearFilter
} from 'https://cdn.skypack.dev/three@0.135.0/build/three.module.js';
import { OrbitControls } from 'https://cdn.skypack.dev/three@0.135.0/examples/jsm/controls/OrbitControls.js';
import {OBJLoader} from 'https://cdn.skypack.dev/three@0.135.0/examples/jsm/loaders/OBJLoader.js';

//0.126.1 0.137.5
import { RoomEnvironment } from 'https://unpkg.com/three@0.126.1/examples/jsm/environments/RoomEnvironment.js';

//--------------------SIMULATION ELEMENTS---------------------

//SIMULATION PARAMETERS GO HERE

function elapsedMs()
{
	return (new Date()).getTime();
}


class Objsim{
	constructor(_frames){
		this.frames = _frames;
		this.currentFrame = 0;
		this.shape = "I HAVE NO SHAPE YET";
	}

	hasLoaded()
	{
		return true;
	}

	nextFrame(){
		var currentFrameTime =  this.frames[this.currentFrame][16]
		while( ( 1000 * currentFrameTime + delay < var_multiplier * (elapsedMs() - startMs)) && (this.currentFrame < this.frames.length-1) )
		{
			this.currentFrame = this.currentFrame + 1;
			currentFrameTime =  this.frames[this.currentFrame][16]
		}


		this.shape.matrix.set( this.frames[this.currentFrame][ 0],this.frames[this.currentFrame][ 1],this.frames[this.currentFrame][ 2],this.frames[this.currentFrame][ 3],
							   this.frames[this.currentFrame][ 4],this.frames[this.currentFrame][ 5],this.frames[this.currentFrame][ 6],this.frames[this.currentFrame][ 7],
							   this.frames[this.currentFrame][ 8],this.frames[this.currentFrame][ 9],this.frames[this.currentFrame][10],this.frames[this.currentFrame][11],
							   this.frames[this.currentFrame][12],this.frames[this.currentFrame][13],this.frames[this.currentFrame][14],this.frames[this.currentFrame][15]);


	}
}

class Box extends Objsim{
	constructor(_name, _width, _height, _depth, _frames, _material){
		super(_frames);
		this.width = _width;
		this.height = _height;
		this.depth = _depth;
		const geometry = new BoxBufferGeometry( this.width, this.height, this.depth);
		const cube = new Mesh( geometry, _material );
		cube.name = _name;
		cube.matrixAutoUpdate = false;
		this.shape = cube;
	}
}

class Ball extends Objsim{
	constructor(_name, _radius, _frames, _material){
		super(_frames);
		this.radius = _radius;
		const geometry = new SphereBufferGeometry( this.radius, 64, 32);
		const sphere = new Mesh( geometry, _material);
		sphere.name = _name;
		sphere.matrixAutoUpdate = false;
		this.shape = sphere;
	}
}

class Cylinder extends Objsim{
	constructor(_name, _radius, _height, _frames, _material){
		super(_frames);
		this.radius = _radius;
		this.height = _height;
		const geometry = new CylinderBufferGeometry( this.radius, this.radius, this.height, 20 );
		const cylinder = new Mesh( geometry, _material);
		cylinder.name = _name;
		cylinder.matrixAutoUpdate = false;
		this.shape = cylinder;
	}
}

class Robot extends Objsim{


	constructor(_objBase, _link, _frames){
		super(_frames);

		this.loadedObjs = 0
		this.totalNumberOfObjects = _link.length+1
		this.links = _link;

		//Function that creates a generic robot
		const base = new Group();
		base.name = "base";

		//Create link groups
		let linkGroupList = [];

		const linkGroup = new Group();
		const axesHelper = new AxesHelper(0.2);

		axesHelper.matrixAutoUpdate = false;
		linkGroup.add(axesHelper)
		linkGroup.name = "link0";

		linkGroupList.push(linkGroup);

		for (let i = 0; i < this.links.length; i++) {

			const linkGroup = new Group();
			const axesHelper = new AxesHelper(0.2);

			axesHelper.matrixAutoUpdate = false;
			linkGroup.add(axesHelper)
			linkGroup.name = "link" + (i+1).toString();

			linkGroup.rotateZ(this.links[i].theta);
			linkGroup.translateZ(this.links[i].d);
			linkGroup.rotateX(this.links[i].alpha);
			linkGroup.translateX(this.links[i].a);

			linkGroup.updateMatrix();

			linkGroupList.push(linkGroup);
		}

		base.add(linkGroupList[0])

		for (let i = 0; i < this.links.length; i++) {
			linkGroupList[i].add(linkGroupList[i + 1])
		}


		//Load 3D models

		this.delta_config = []

		this.loadObj(_objBase,"base", "base", false)
		for (let i = 0; i < this.links.length-1; i++) {
			this.loadObj(this.links[i].model3d,"link" + (i+1).toString(), "link" + (i).toString(), false)
		}
		this.loadObj(this.links[this.links.length-1].model3d,"link" + (this.links.length).toString(), "link" + (this.links.length).toString(), true)

		this.shape = base;
		this.shape.matrixAutoUpdate = false;

		for (let i = 0; i < this.links.length; i++) {
			if(this.links[i].jointType == 0)
			{
				this.delta_config.push(this.shape.getObjectByName("link" + (i).toString()).rotation.z)
			}
			if(this.links[i].jointType == 1)
			{
				this.delta_config.push(this.shape.getObjectByName("link" + (i).toString()).position.z)
			}
		}
	}

	loadObj(obj, name1, name2, axisVisible) {
		const manager = new LoadingManager();
		const objLoader = new OBJLoader(manager);
		objLoader.load(obj.url, (root) => {
			root.scale.set(obj.scale, obj.scale, obj.scale);
			root.applyMatrix4(obj.matrix)

			root.traverse(function (child) {
				if (child instanceof Mesh) {
					child.material = obj.mesh_material;
				}
			});


			this.shape.getObjectByName(name1).getObjectByProperty("type", "AxesHelper").visible = axisVisible;
			this.shape.getObjectByName(name2).add(root);
			this.loadedObjs += 1
		});
	}

	hasLoaded()
	{
		return this.loadedObjs == this.totalNumberOfObjects;
	}

	//Function that updates frames
	nextFrame(){

		var currentFrameTime =  this.frames[this.currentFrame][17]
		while( ( 1000 * currentFrameTime + delay < var_multiplier * (elapsedMs() - startMs)) && (this.currentFrame < this.frames.length-1) )
		{
			this.currentFrame = this.currentFrame + 1;
			currentFrameTime =  this.frames[this.currentFrame][17]
		}


		//setting robot position
		this.shape.matrix.set( this.frames[this.currentFrame][ 0],this.frames[this.currentFrame][ 1],this.frames[this.currentFrame][ 2],this.frames[this.currentFrame][ 3],
							   this.frames[this.currentFrame][ 4],this.frames[this.currentFrame][ 5],this.frames[this.currentFrame][ 6],this.frames[this.currentFrame][ 7],
							   this.frames[this.currentFrame][ 8],this.frames[this.currentFrame][ 9],this.frames[this.currentFrame][10],this.frames[this.currentFrame][11],
							   this.frames[this.currentFrame][12],this.frames[this.currentFrame][13],this.frames[this.currentFrame][14],this.frames[this.currentFrame][15]);
		//setting robot configuration

		if(this.frames[this.currentFrame][16] != undefined){
			let linkName = "";

			for(let i = 0; i < this.links.length; i++){
				linkName = "link" + (i).toString();

				if(this.links[i].jointType == 0){
					this.shape.getObjectByName(linkName).rotation.z = this.frames[this.currentFrame][16][i] + this.delta_config[i];
				}else if (this.links[i].jointType == 1){
					this.shape.getObjectByName(linkName).position.z = this.frames[this.currentFrame][16][i] + this.delta_config[i];
				}
			}
		}


	}
}

//------------------------------------------------------------

//--------------- BASIC ELEMENTS OF ANY SCENE ----------------
Object3D.DefaultUp = new Vector3(0,0,1); //Pointing Z axis up
const canvas = document.querySelector('#scene');// Selecting canvas

const scene = new Scene();//Instantiate the Scene

scene.background = new Color('white');//Set background color
const camera = new PerspectiveCamera(35, canvas.clientWidth/canvas.clientHeight, 0.1, 100);//Instantiate a camera
camera.position.set(4, 4, 3);//Put camera in its place
var ambientLight = new HemisphereLight('white','darkslategrey', 3);//Instantiate Ambient light
const controls = new OrbitControls(camera, canvas);	//Instantiate orbit controls
controls.target.set(0, 0, 0);//Point camera at the origin
const renderer = new WebGLRenderer({canvas, antialias: true});//Instantiate renderer
renderer.physicallyCorrectLights = true;//Enable physically Correct Lights
renderer.setSize(canvas.clientWidth, canvas.clientHeight);//Set render size
renderer.setPixelRatio(window.devicePixelRatio);//Set pixel ratio



let sceneElements = [];
const axesHelper = new AxesHelper( 0.5 ); //Create axis helper
axesHelper.renderOrder = 1;
const gridHelper = new GridHelper( 3, 6);//Create grid helper
gridHelper.rotation.x = 3.14/2;

scene.add(ambientLight);
scene.add( axesHelper );
scene.add( gridHelper );
//------------------------------------------------------------

//--------------- ADDING ELEMENTS TO THIS SCENE ---------------
//USER INPUT GOES HERE

// add stuff to the scene
for(let i = 0; i < sceneElements.length; i++){
	scene.add(sceneElements[i].shape);
}


//------------------------------------------------------------

//Wait until all objects loaded

var loadedAll = false;


//Get starting ms
var startMs =  elapsedMs()


//-------------------- THE ANIMATION LOOP -------------------
renderer.setAnimationLoop(() => {


	loadedAll = true;

	for(let i = 0; i < sceneElements.length; i++){
		loadedAll = loadedAll && sceneElements[i].hasLoaded();
	}

	if(!loadedAll)
	{
		startMs =  elapsedMs()
	}



	//MAGIC HAPPENS HERE!!!
	for(let i = 0; i < sceneElements.length; i++){
		sceneElements[i].nextFrame();
	}



	renderer.render(scene, camera);
});
//------------------------------------------------------------
