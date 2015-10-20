/* 
 * The MIT License
 *
 * Copyright 2015 yl9.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

function relative_angle(observerX, observerY, observeeX, observeeY) {

    var rx = observerX;
    var ry = observerY;
    var ex = observeeX;
    var ey = observeeY;
    var dx = rx - ex;
    var dy = ry - ey;

    /*
     if (Math.abs(dx) > poolWidth / 2) {
     if (rx < ex) {
     dx = (rx + poolWidth) - ex;
     } else if (rx > ex) {
     dx = rx - (ex + poolWidth);
     }
     }
     if (Math.abs(dy) > poolHeight / 2) {
     if (ry < ey) {
     dy = (ry + poolHeight) - ey;
     } else if (ry > ey) {
     dy = ry - (ey + poolHeight);
     }
     }
     */


    //var deg = Math.atan2((-observerY + observeeY), (-observerX + observeeX));
    var deg = Math.atan2(-dy, -dx);
    if (deg < 0) {
        deg += Math.PI * 2;
    }
    return deg;
}

function getVisibleSize() {
    var myWidth = 0, myHeight = 0;
    if (typeof (window.innerWidth) === 'number') {
        //Non-IE
        myWidth = window.innerWidth;
        myHeight = window.innerHeight;
    } else if (document.documentElement && (document.documentElement.clientWidth || document.documentElement.clientHeight)) {
        //IE 6+ in 'standards compliant mode'
        myWidth = document.documentElement.clientWidth;
        myHeight = document.documentElement.clientHeight;
    } else if (document.body && (document.body.clientWidth || document.body.clientHeight)) {
        //IE 4 compatible
        myWidth = document.body.clientWidth;
        myHeight = document.body.clientHeight;
    }
    var result = {width: myWidth, height: myHeight};
    return result;
}
var winSize = getVisibleSize();
var winWidth = 1024;//winSize.width;
var winHeight = 768;//winSize.height;

var poolWidth = 1024;
var poolHeight = 768;

var cellpopulation = 50;
var foodamount = 50;

var hiddenLayerNum = 1;
var cellCount = 0;
var foodCount = 0;

var minRadius = 3;
var maxRadius = 6;

var maxSpeed = 3.5;
var minSpeed = 0;

var cellgrid;
var foodgrid;


var reposition = function (cell) {
    if (cell.stamina <= 0) {
        return;
    }
    var speed = cell.speed;
    cell.stamina -= Math.pow(speed * 2, 1.3) * 0.01;
    var newx = cell.x + Math.cos(cell.r) * cell.speed;
    var newy = cell.y + Math.sin(cell.r) * cell.speed;
    var on_h_edge = false;
    var on_v_edge = false;
    if (newx >= poolWidth - 1) {
        newx = poolWidth - 1;
        on_v_edge = true;

    }
    if (newx <= 1) {
        newx = 1;
        on_v_edge = true;
    }

    if (newy >= poolHeight - 1) {
        newy = poolHeight - 1;
        on_h_edge = true;
    }
    if (newy <= 1) {
        newy = 1;
        on_h_edge = true;
    }

    if (on_h_edge === true) {
        cell.y = newy;
        cell.speed = 0;
        return;
    } else if (on_v_edge === true) {
        cell.x = newx;
        cell.speed = 0;
        return;
    } else {
        cell.x = newx;
        cell.y = newy;
    }


};

var accelerate = function (cell) {
    cell.speed += 0.2;
    if (cell.speed >= maxSpeed) {
        cell.speed = maxSpeed;
    }
};

var decelerate = function (cell) {
    cell.speed -= 0.2;
    if (cell.speed <= minSpeed) {
        cell.speed = minSpeed;
    }
};

var turnLeft = function (cell) {
    cell.r = (cell.r - 2 * Math.PI / 12 / (cell.speed + 0.01)) % (Math.PI * 2); //turn left for 360/24 degree

};

var turnRight = function (cell) {
    cell.r = (cell.r + 2 * Math.PI / 12 / (cell.speed + 0.01)) % (Math.PI * 2); //turn left for 360/24 degree

};

var defend = function (cell) {
    cell.is_defending = true;
};

function getNeighborObjs(obj, grid) {
    var result = [];
    var pos = grid.xyTocr(obj.x, obj.y);
    var col = pos[0];
    var row = pos[1];
    var grids = grid.getNeighborGrids(col, row);
    //console.log("test:",grids);
    for (var i = 0; i < grids.length; i++) {
        //for (var j = 0; j < grids[i].length; j++) {
        var neighbor = grids[i];
        result = result.concat(neighbor);
        //}
    }
    //if (framecount % 60 === 0)
    //console.log(grids.length, col, row, result);
    return result;
}
var distance = function (x1, y1, x2, y2) {
    //var dx = Math.min(Math.abs(x2 - x1), Math.abs(x2 + poolWidth - x1), Math.abs(x2 - (x1 + poolWidth)));
    //var dy = Math.min(Math.abs(y2 - y1), Math.abs(y2 + poolHeight - y1), Math.abs(y2 - (y1 + poolHeight)));
    var dx = x1 - x2;
    var dy = y1 - y2;

    return Math.sqrt(dx * dx + dy * dy);
};

var attack = function (cell) {


    cell.stamina -= 1;
    var neignbors = getNeighborObjs(cell, cellgrid);


    for (var i = 0; i < neignbors.length; i++) {
        var vcell = neignbors[i];
        var dist = distance(cell.x, cell.y, vcell.x, vcell.y);
        if (dist < cell.radius + vcell.radius) {
            if (cell.id !== vcell.id && vcell.stamina > 0) {
                //console.log("cell id=" + cell.id + " attacked cell id=" + vcell.id);
                vcell.is_attacked = true;
                vcell.attackers.push(cell);
            }

        }
    }

};

var cellDefault = {
    inputNum: 7,
    outputNum: 6,
    //list of action functoins
    //3 move function allows the cell to move at 3 different speed
    actions: [
        accelerate,
        decelerate,
        turnLeft,
        turnRight,
        attack,
        defend

                //defend,
                // attack
    ]


};

function Food(x, y, energy) {
    this.id = foodCount++;

    // console.log("new food created id=" + this.id);
    if (typeof x === "undefined" || typeof y === "undefined") {
        this.x = Math.random() * poolWidth;
        this.y = Math.random() * poolHeight;
    } else {
        this.x = x;
        this.y = y;
    }
    //for grid compatibility
    this.energy = 80;// 100.0 * Math.random();
    if (typeof energy !== "undefined") {
        this.energy = energy;
    }

    this.radius = minRadius;
    this.color = [
        Math.floor(100 / 100 * this.energy),
        Math.floor(255 / 100 * this.energy),
        Math.floor(100 / 100 * this.energy)];
    this.lastx = this.x;
    this.lasty = this.y;

    this.is_eaten = false;
}

function fixColorVec(color) {
    if (color < 0) {
        return 0;
    } else if (color > 255) {
        return 255;
    } else {
        return Math.floor(color);
    }
}

function Cell(genome, x, y, r) {
    // x , y coordinate
    // r rotation
    //neuro , the neuroo network of this cell

    //id is used to identify cell.
    this.id = cellCount++;
    this.color = [255, 255, 255];
    this.radius = 8;
    this.sensor_dist = 100;
    if (typeof x === "undefined" || typeof y === "undefined" || typeof r === "undefined") {
        this.x = Math.random() * poolWidth;
        this.y = Math.random() * poolHeight;
        this.r = Math.random() * Math.PI;
    } else {
        this.x = x;
        this.y = y;
        this.r = r;
    }
    //coordinates of this cell in last frame
    this.lastx = this.x;
    this.lasty = this.y;
    this.speed = 0;

    //this.is_moving = false;

    this.is_attacked = false;
    this.attackers = [];
    this.is_defending = false;

    //this.is_running = false;

    this.is_dead = false;

    this.is_eating = false;
    this.food_in_mouth = [];

    if (typeof genome === "undefined") {
        this.genome = null;
    } else {
        this.genome = genome;
        this.genome.fitness = 0;
    }
    this.stamina = 800;
    this.life = 100;
    this.fitness = 0;


}
Cell.prototype.init = function () {
    this.stamina = 800;
    this.life = 100;
    this.fitness = 0;
    if (typeof this.genome !== "undefined")
        this.genome.fitness = 0;

    this.is_moving = false;

    this.is_attacked = false;
    this.attackers = [];
    this.is_defending = false;
    this.is_running = false;

    this.is_dead = false;

    this.is_eating = false;
    this.food_in_mouth = [];
};

function getNeighborObjFromGrid(cell, grid) {

    var objAngle = 0;
    var cr = 0, cg = 0, cb = 0;

    var neighbors = getNeighborObjs(cell, grid);
    var objDist = cell.sensor_dist + cell.radius;
    var minObjDist = objDist;
    var nearest_obj;
    if (cell.id === debugid) {
        // if (framecount % 30 === 0)
        //console.log(neighbors);

    }
    for (var i = 0; i < neighbors.length; i++) {
        var dist = distance(cell.x, cell.y, neighbors[i].x, neighbors[i].y);
        if (dist < minObjDist) {
            if (typeof neighbors[i].is_dead !== "undefined") {

                if (neighbors[i].is_dead === true) {
                    continue;
                }
                if (neighbors[i].id === cell.id) {
                    continue;
                }
            }
            if (typeof neighbors[i].is_eaten !== "undefined") {
                if (neighbors[i].is_eaten === true) {
                    continue;
                }
            }
            minObjDist = dist;
            nearest_obj = neighbors[i];
        }
    }
    if (typeof nearest_obj !== "undefined") {
        if (minObjDist < cell.radius + cell.sensor_dist) {
            objDist = minObjDist;
            objAngle = relative_angle(cell.x, cell.y, nearest_obj.x, nearest_obj.y) - cell.r;

            if (objAngle < -Math.PI) {
                objAngle += Math.PI * 2;
            } else if (objAngle > Math.PI) {
                objAngle -= Math.PI * 2;
            }

            cr = nearest_obj.color[0];
            cg = nearest_obj.color[1];
            cb = nearest_obj.color[2];

        }
    }
    var result = {obj: nearest_obj, data: [objDist, objAngle, cr, cg, cb]};
    //if (cell.id === 0) {
    //     if (framecount % 60 === 0)
    //       console.log(neighbors, result);
    //}

    return result;
}

Cell.prototype.sense = function () {
    var wallDist = Math.min(this.x, this.y, poolWidth - this.x, poolHeight - this.y);
    var wallSensor;
    if (wallDist <= 2) {
        wallSensor = 1;
    } else {
        wallSensor = 0;
    }
    var v1 = getNeighborObjFromGrid(this, cellgrid);
    var v2 = getNeighborObjFromGrid(this, foodgrid);
    var rv1 = [v1.data[0], v1.data[1]];
    var rv2 = [v2.data[0], v2.data[1]];

    var foodDist = v2.data[0];
    var food = v2.obj;
    if (foodDist < this.radius) {
        if (food.is_eaten === false) {
            this.is_eating = true;
            food.is_eaten = true;
            this.food_in_mouth = [food];
        }
    }



    var result = rv1.concat(rv2).concat([this.speed, this.stamina, wallSensor]);//.concat(wallDist);
    if (this.id === debugid) {
        if (this.is_dead == false) {
            //if (framecount % 30 === 0)
            //console.log(result);
        }
    }
    // if (this.id === 0 && this.is_dead === false) {
    //      if (framecount % 15 === 0) {
    //          console.log("" + result);
    //     }
    //  }

    return result;
};
Cell.prototype.toVector = function () {
    var result = [this.color[0], this.color[1]];
    result = result.concat(this.neuro.toVector());
    return result;
};



function contains(list, item) {
    //console.log(list);
    for (var i = 0; i < list.length; i++) {
        if (list[i] === item) {
            return i;
        }
    }
    //console.log("return -1");
    return -1;
}

function Gridmap(col, row) {
    this.colnum = col;
    this.rownum = row;
    this.gridwidth = poolWidth / col;
    this.gridheight = poolHeight / row;
    this.gridmap = new Array(col);

    //create grdimap [x][y]
    for (var x = 0; x < col; x++) {
        this.gridmap[x] = new Array(row);
        for (var y = 0; y < row; y++) {

            //console.log(x + "," + y);
            this.gridmap[x][y] = [];

        }

    }

}
Gridmap.prototype.xyTocr = function (x, y) {
    var c = Math.floor(x / this.gridwidth);
    var r = Math.floor(y / this.gridheight);
    return [c, r];
};
Gridmap.prototype.updateObjPos = function (obj) {
    //console.log(obj);
    var pos = this.xyTocr(obj.x, obj.y);
    var c = pos[0];
    var r = pos[1];
    //console.log(cell.x + "," + cell.y);
    //console.log(this.gridwidth + "," + this.gridheight);
    //console.log(c + "," + r);
    //console.log(obj);

    //console.log("index=" + index);
    var last_pos = this.xyTocr(obj.lastx, obj.lasty);
    var lastc = last_pos[0];
    var lastr = last_pos[1];
    var last_index = contains(this.gridmap[lastc][lastr], obj);

    if (last_index === -1) {
        //if last and current frame don't have the cell, it means the cell was just created.
        //simply add cell
    } else {
        //if the movement of this cell has changed the grid which the cell should belong
        this.gridmap[lastc][lastr].splice(last_index);

    }

    var index = contains(this.gridmap[c][r], obj);
    if (index === -1) {
        this.gridmap[c][r].push(obj);
        // console.log("pushed");
    } else {

    }
};
Gridmap.prototype.getObjs = function (col, row) {
    return this.gridmap[col][row];
};
Gridmap.prototype.getNeighborGrids = function (col, row) {
    var result = [];
    var n = 3;
    var sc = col - n;
    var ec = col + n;
    var sr = row - n;
    var er = row + n;
    if (sc <= 0)
        sc = 0
    if (ec >= this.colnum)
        ec = this.colnum - 1

    if (sr <= 0)
        sr = 0
    if (er >= this.rownum)
        er = this.rownum - 1


    for (var c = sc; c <= ec; c++) {
        for (var r = sr; r <= er; r++) {
            var real_c = c;// (c + this.colnum) % this.colnum;
            var real_r = r;// (r + this.rownum) % this.rownum;
            result.push(this.gridmap[real_c][real_r]);
        }
    }
    return result;
};
Gridmap.prototype.removeObj = function (obj) {
    var lastPos = this.xyTocr(obj.lastx, obj.lasty);
    var lastc = lastPos[0];
    var lastr = lastPos[1];
    var index = contains(this.gridmap[lastc][lastr], obj);
    if (index !== -1) {
        this.gridmap[lastc][lastr].splice(index, 1);
    }

    var Pos = this.xyTocr(obj.lastx, obj.lasty);
    var c = Pos[0];
    var r = Pos[1];
    var index = contains(this.gridmap[c][r], obj);
    if (index !== -1) {
        this.gridmap[c][r].splice(index, 1);
    }


};

function createCells(num) {
    var result = [];
    for (var i = 0; i < num; i++) {
        result.push(new Cell());
    }
    return result;
}

function createFoods(num) {
    var result = [];
    for (var i = 0; i < num; i++) {
        result.push(new Food());
    }
    return result;
}

var cellgrid = new Gridmap(28, 16);
var foodgrid = new Gridmap(28, 16);
var foods = [];
function initFoods() {
    // for (var i = 0; i < foods.length; i++) {
    //   foodgrid.removeObj(foods[i]);
    // }
    foods = createFoods(foodamount);
    for (var i = 0; i < foodamount; i++) {
        foodgrid.updateObjPos(foods[i]);

    }
}