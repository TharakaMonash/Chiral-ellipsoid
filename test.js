let WIDTH = 5;
 let LENGTH = 3;
 let topAndBottom = new Array(0)
 let middles = new Array(0);
 for( let i = 0; i < WIDTH; i++){
 topAndBottom.push("*");
 }
 let topAndBottomStr = topAndBottom.join('');
 console.log(topAndBottomStr);

 middles.push("*");
 for( let k = 0; k < (WIDTH - 2);k++){
 middles.push(" ");
 }
 middles.push("*");

 middlesStr = middles.join('');
 for( let j = 0; j< (LENGTH - 2); j++){
 console.log(middlesStr);
 }

 console.log(topAndBottomStr)
