/*
 Wave sliders

*/


let line;
let sineGen;
let envs;
let showGraphs;
let dur;
let counter;

let ranges;

let ampSlider;
let incSlider;
let freqSlider;
let warpSlider;
let warpCentroidSlider;


function setup()
{




  createCanvas(512,512);
  frameRate(60);

  let vertices = [];
  for(let i = 0; i < 100; i++)
  {
    vertices.push(createVector(-width/2 + i * width/100,0));
  }

  sineGens = [];

  line = new SimpleLine(vertices);
  sineGen = new SineGen(1,100,0);



  envs = {};
  ranges = {};
  //amp

  ampSlider = createSlider(0, height/2).position(30,20);
  incSlider = createSlider(PI/500, PI/10, PI/500, PI/500).position(30,40);
  freqSlider = createSlider(0.001,5,0.001,0.001).position(30,60);
  warpSlider = createSlider(0.01,1,1,0.01).position(30,80);
  warpCentroidSlider = createSlider(0,1,0.5,0.01).position(30,100);


  dur = 20;
  counter = 0;




}

function draw()
{
  background(255);
  noFill();

  text("amp", 10,20);
  text("inc", 10,40);
  text("freq", 10,60);
  text("warp " + sineGen.warp , 10,80);
  text("warpCentroid" + sineGen.warpCentroid , 10,100);

  sineGen.amp =  ampSlider.value(); //map( envs.amp.lin_value(p), 0,1,ranges['amp'].min, ranges['amp'].max);
  sineGen.freq =  freqSlider.value();
  let inc = incSlider.value();
  sineGen.warp = pow(warpSlider.value(),2) * 3;
  sineGen.warpCentroid = warpCentroidSlider.value();
  sineGen.update(inc);




  //draw the shape
  stroke(0);
  translate(width/2,height/2);

  beginShape();
  for(let i = 0; i < 512; i++)
  {
    let v = line.calcVertex(i/512);
    let y = sineGen.value(i/512);
    vertex(v.x,v.y + y);
  }
  endShape();





}

function keyPressed()
{
  if(key == ' ')
  {
    counter = millis();
  }
  else if(key =='v')
  {
    showGraphs = !showGraphs;
  }
}
