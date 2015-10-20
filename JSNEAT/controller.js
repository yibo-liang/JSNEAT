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

var cell_alive_count = cellpopulation;
var lastGen = [];

var Info = function (gen, foodnum, cellnum) {
    this.gen = gen;
    this.foodnum = foodnum;
    this.cellnum = cellnum;
};

function cellAct(cellData) {
    for (var i = 0; i < cellData.length; i++) {
        var cell = cellData[i];
        var input = cell.sense();
        generateNeuroNetwork(cell.genome);

        var output = evaluateNeuroNetwork(cell.genome.network, input);
        if (cellData[i].id === debugid && framecount % 30 === 0 && !is_fastmode) {
            //console.log(input, output);
        }

        //continue;

        for (var j = 0; j < output.length; j++) {
            var f = cellDefault.actions[j];
            //console.log("cell id=" + cell.id + " in action " + f.name);
            if (output[j] === true && cell.stamina > 0)
                f(cell);
        }
    }
}

function cellUpdate(cellData) {
    //console.log("cell update len="+cellData.length);
    var liveCells = [];
    var deadCells = [];
    for (var i = 0; i < cellData.length; i++) {
        var cell = cellData[i];
        //console.log(cell);
        if (true) {
            cell.lastx = cell.x;
            cell.lasty = cell.y;
            reposition(cell);
            //cell.is_moving = false;

        }
        if (cell.is_attacked) {
            for (var a = 0; a < cell.attackers.length; a++) {
                if (cell.is_defending) {
                    cell.life -= 5;
                } else {
                    cell.life -= 10;
                }
                //console.log("cell id=" + cell.id + " is attacked by id=" + cell.attackers[a].id);
                if (cell.life <= 0) {
                    // console.log("cell id=" + cell.id + " is killed by id=" + cell.attackers[a].id);
                    cell.attackers[a].color[0] *= 0.9;
                    break;
                }

            }
            cell.is_attacked = false;
            cell.is_defending = false;
            cell.attackers = [];
        }
        if (cell.is_eating) {
            cell.is_eating = false;
            var energy = cell.food_in_mouth[0].energy;
            cell.food_in_mouth = [];
            cell.stamina += energy;
            cell.life += energy * 0.8;
            if (cell.life > 100)
                cell.life = 100;
        }

        cell.life -= 0.08; // each frame the cell loses 0.1 life, out of 100
        if (cell.stamina <= 0) {
            //cell.life -= 0.12;
        }
        if (cell.life <= 0) {
            cell.is_dead = true;
            cell_alive_count -= 1;
            //console.log("cell id=" + cell.id + " ran out of life, living cell = " + cell_alive_count);
            cellgrid.removeObj(cell);
            deadCells.push(cell);
        } else {
            //each frame the cell is alive, its fitness increases
            var hp = cell.life / 100;
            if (hp > 1)
                hp = 1;
            cell.color[2] = fixColorVec(hp * 255);
            //cell.radius = Math.floor(minRadius + (maxRadius - minRadius) * hp);
            cellgrid.updateObjPos(cell);
            cell.fitness++;
            cell.genome.fitness = cell.fitness;
            liveCells.push(cell);
        }

    }
    return {lc: liveCells, dc: deadCells};
}

var foodUpdate = function (foodData, deadCells) {
    var result = [];
    var is_dfood = true;
    if (is_dfood === true) {
        for (var j = 0; j < deadCells.length; j++) {
            // when cell dies, it turns to food with more energy in it.
            var newfood = new Food(deadCells[j].x, deadCells[j].y, 85);
            newfood.color = [255, 0, 0];
            result.push(newfood);
            //console.log("dead cell id=" + deadCells[j].id + " becomes food, living cell = " + cell_alive_count);
        }
    }

    for (var i = 0; i < foodData.length; i++) {
        if (foodData[i].is_eaten === false) {
            result.push(foodData[i]);

        } else {
            foodgrid.removeObj(foodData[i]);
        }
    }

    for (var i = 0; i < result.length; i++) {
        foodgrid.updateObjPos(result[i]);
    }

    return result;
};

var is_fastmode = false;
var fastmode_gen = 5;
var generation = 0;

function doFrame() {
    var lastcells = cellUpdate(cells);
    cells = lastcells.lc;
    var deadCells = lastcells.dc;
    lastGen = lastGen.concat(deadCells);


    foods = foodUpdate(foods, deadCells);
    if (framecount % 60 == 0) {
        // console.log(cells[0].neuro);
    }

    if (cell_alive_count === 0) {
        cellgrid = new Gridmap(28, 16);
        foodgrid = new Gridmap(28, 16);
        //console.log(lastGen);

        newGeneration();

        console.log("Top Fitness=" + pool.maxFitness, "staleness=" + pool.staleness);

        cells = [];
        console.log(pool);
        for (var i = 0; i < pool.species.length; i++) {
            for (var j = 0; j < pool.species[i].genomes.length; j++) {
                var nc = new Cell(pool.species[i].genomes[j]);
                cells.push(nc);
            }
        }
        //Population = cells.length;
        //cellpopulation = Population;

        lastGen = [];
        cell_alive_count = cells.length;

        console.log("new gen, new cell amount=" + cells.length);
        generation++;
        initFoods();
        fastmode_gen--;
        //if (fastmode_gen === 0)
        //is_fastmode = false;

    }
    info.cellnum = cells.length;
    info.foodnum = foods.length;
    info.gen = generation;
    d3.select("body").select(".panel").selectAll("p").data([info]).text(function (d) {
        return  "gen:" + d.gen
                + ", food n=" + d.foodnum
                + ", cell n=" + d.cellnum;
    });

    if (cells.length === 1) {
        if (framecount % 60 === 0) {
            //console.log(cells[0].sense());
        }
    }
    cellAct(cells);
    var cellDataBinding = dataContainer.selectAll("container.cell").data(cells);
    var foodDataBinding = dataContainer.selectAll("container.food").data(foods);

    cellDataBinding.enter().append("container").classed("cell", true);
    foodDataBinding.enter().append("container").classed("food", true);

    function f() {
        if (!is_fastmode)
            render(cellDataBinding, foodDataBinding);
    }
    f();
}

function doFrameFast() {
    if (is_fastmode)
        setTimeout(function () {
            doFrame();
            doFrameFast();
        }, 0);
}

function speedToggle() {
    is_fastmode = !is_fastmode;
    if (is_fastmode)
        fastmode_gen = 20;
}

var stepframe = function () {

    if (is_fastmode) {
        doFrameFast();
        return;
    } else {
        doFrame();
    }




};