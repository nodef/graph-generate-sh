const fs = require('fs');
const os = require('os');
const path = require('path');

const RGRAPH = /^(?:Loading|Using) graph .*\/(.*?)\.\w+ \.\.\./m;
const RGRAPD = /^order: (\d+) size: (\d+) {}/m;
const RTEMPE = /^- temporal-edges: (\d+)/m;
const RSPACE = /^---|^$/m;
const RFIELD = /^- (.+?):\s+(\S+)/m;




// *-ARRAY
// -------

function sumArray(x) {
  var a = 0;
  for (var i=0, I=x.length; i<I; i++)
    a += x[i];
  return a;
}

function avgArray(x) {
  return x.length? sumArray(x)/x.length : 0;
}




// *-FILE
// ------

function readFile(pth) {
  var d = fs.readFileSync(pth, 'utf8');
  return d.replace(/\r?\n/g, '\n');
}

function writeFile(pth, d) {
  d = d.replace(/\r?\n/g, os.EOL);
  fs.writeFileSync(pth, d);
}




// *-CSV
// -----

function writeCsv(pth, rows) {
  var cols = Object.keys(rows[0]);
  var a = cols.join()+'\n';
  for (var r of rows)
    a += [...Object.values(r)].map(v => `"${v}"`).join()+'\n';
  writeFile(pth, a);
}




// *-LOG
// -----

function readLogLine(ln, data, state) {
  if (RGRAPH.test(ln)) {
    var [, graph] = RGRAPH.exec(ln);
    if (!data.has(graph)) data.set(graph, []);
    state = {graph, temporal_edges: 0, fields: {}};
  }
  else if (RTEMPE.test(ln)) {
    var [, temporal_edges] = RTEMPE.exec(ln);
    state.temporal_edges = parseFloat(temporal_edges);
  }
  else if (RSPACE.test(ln)) {
    if (Object.keys(state.fields).length>0) {
      data.get(state.graph).push(Object.assign({
        graph:          state.graph,
        temporal_edges: state.temporal_edges
      }, state.fields));
    }
    state.fields = {};
  }
  else if (RFIELD.test(ln)) {
    var [, name, value] = RFIELD.exec(ln);
    name = name.replace(/-/g, '_');
    state.fields[name] = value;
  }
  return state;
}

function readLog(pth) {
  var text = readFile(pth);
  var lines = text.split('\n');
  var data = new Map();
  var state = null;
  for (var ln of lines)
    state = readLogLine(ln, data, state);
  return data;
}




// PROCESS-*
// ---------

function processCsv(data) {
  var a = [];
  for (var rows of data.values())
    a.push(...rows);
  return a;
}




// MAIN
// ----

function main(cmd, log, out) {
  var data = readLog(log);
  if (path.extname(out)==='') cmd += '-dir';
  switch (cmd) {
    case 'csv':
      var rows = processCsv(data);
      writeCsv(out, rows);
      break;
    case 'csv-dir':
      for (var [graph, rows] of data)
        writeCsv(path.join(out, graph+'.csv'), rows);
      break;
    default:
      console.error(`error: "${cmd}"?`);
      break;
  }
}
main(...process.argv.slice(2));
