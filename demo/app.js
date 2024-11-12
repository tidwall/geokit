function run(opts) {
////////////////////////////////////////////////////////////////////////////////

opts = opts || {};
opts.name = opts.name || 'undefined';
opts.links = opts.links || [ 
    'direct', 
    'inverse', 
    'track',
    'intercept', 
    'intersect', 
    'nearest',
    'equibox',
    'equipath',
] ;
opts.show = opts.show || [ 'a', 'b', 'c', 'd' ];
opts.draw = opts.draw || function(){};
opts.about = opts.about || '';

////////////////////////////////////////////////////////////////////////////////

// write the entire
document.body.innerHTML += `
<div id="controls"></div>
<div id="map"></div>
<div id="topright-links">
    <a href="https://github.com/tidwall/geokit" id="github-link" style="filter: invert(0); display: inline-block; width: 32px; height: 32px ">
        <svg width="32" height="32" aria-hidden="true" viewBox="0 0 16 16" version="1.1" data-view-component="true">
            <path d="M8 0c4.42 0 8 3.58 8 8a8.013 8.013 0 0 1-5.45 7.59c-.4.08-.55-.17-.55-.38 0-.27.01-1.13.01-2.2
                0-.75-.25-1.23-.54-1.48 1.78-.2 3.65-.88 3.65-3.95 0-.88-.31-1.59-.82-2.15.08-.2.36-1.02-.08-2.12
                0 0-.67-.22-2.2.82-.64-.18-1.32-.27-2-.27-.68 0-1.36.09-2 .27-1.53-1.03-2.2-.82-2.2-.82-.44
                1.1-.16 1.92-.08 2.12-.51.56-.82 1.28-.82 2.15 0 3.06 1.86 3.75 3.64 3.95-.23.2-.44.55-.51
                1.07-.46.21-1.61.55-2.33-.66-.15-.24-.6-.83-1.23-.82-.67.01-.27.38.01.53.34.19.73.9.82
                1.13.16.45.68 1.31 2.69.94 0 .67.01 1.3.01 1.49 0 .21-.15.45-.55.38A7.995 7.995 0 0 1 0
                8c0-4.42 3.58-8 8-8Z"></path>
        </svg>
    </a>
</div>
<div id="bottom-about">
    <div align=center id=best-viewed><b>best viewed on desktop</b></div>
    This site is a demonstration of the various geodesic calculations provided by the GeoKit C library.<br>
    To learn more about GeoKit visit <a href="https://github.com/tidwall/geokit">github.com/tidwall/geokit</a>.<br>
    © <span id="jbcopyright">Copyright 2023, &nbsp;</span>Josh Baker
</div>
`

let el;

el = document.getElementById('controls');
if (!el) throw new Error("missing 'controls' element");
el.innerHTML = `
    <div id="name">` + opts.name + `</div>
    <div id="nav"></div>
    <div id="about"></div>
    <div id="info-details">
        <div id="sphere_info">
            <div class="shapehdr hblue"><div id="toggle-sphere"></div><pre>   Sphere</pre></div>
            <pre id="vals_sphere"></pre>
            <pre id="vals_sphere_min"></pre>
        </div>
        <div id="ellipsoid_info">
            <div class="shapehdr hred"><div id="toggle-ellipsoid"></div><pre>   Ellipsoid</pre></div>
            <pre id="vals_ellipsoid"></pre>
            <pre id="vals_ellipsoid_min"></pre>
        </div>
    </div>
`;

el = document.getElementById('nav');
for (let i = 0; i < opts.links.length; i++) {
    if (i > 0) {
        el.innerHTML += '&nbsp;• ';
    }
    el.innerHTML += '<a href="' + opts.links[i] + '.html">' + 
        opts.links[i] + '</a>';
}

let isDisabled = {
    'sphere': sessionStorage.getItem("sphereDisabled") == "true",
    'ellipsoid': sessionStorage.getItem("ellipsoidDisabled") == "true",
}

el = document.getElementById('map');
if (!el) throw new Error("missing 'map' element");

const defPosn = document.body.offsetWidth < 640 ? `[70,0,0]` : `[0,0,1]`;
const stateName = opts.name + ":state";
const stateNamePosn = "shared:state:posn";
let init = JSON.parse(sessionStorage.getItem(stateName) || `
    {
        "isdef": true,
        "posn": `+defPosn+`,
        "a": [19.31114, -98.43750],
        "b": [-28.82259, 40.07813],
        "c": [58.81374, 7.03125],
        "d": [-59.36828, -74.53125]
    }
`);
let sposn = JSON.parse(sessionStorage.getItem(stateNamePosn) || defPosn);
if (sposn) {
    
    init.posn = sposn;
}
if (init.isdef) {
    delete init.isdef;
    if (opts.init) {
        if (opts.init.a) {
            if (opts.init.a.lat) init.a[0] = opts.init.a.lat;
            if (opts.init.a.lon) init.a[1] = opts.init.a.lon;
        }
        if (opts.init.b) {
            if (opts.init.b.lat) init.b[0] = opts.init.b.lat;
            if (opts.init.b.lon) init.b[1] = opts.init.b.lon;
        }
        if (opts.init.c) {
            if (opts.init.c.lat) init.c[0] = opts.init.c.lat;
            if (opts.init.c.lon) init.c[1] = opts.init.c.lon;
        }
        if (opts.init.d) {
            if (opts.init.d.lat) init.d[0] = opts.init.d.lat;
            if (opts.init.d.lon) init.d[1] = opts.init.d.lon;
        }
        if (opts.init.center) {
            if (opts.init.center.lat) init.posn[0] = opts.init.center.lat;
            if (opts.init.center.lon) init.posn[1] = opts.init.center.lon;
        }
        if (opts.init.zoom) init.posn[2] = opts.init.zoom;
    }
}


let autoStoring = false;
const saveState = function(a, b, c, d, map) {
    init = {
        "posn": [ map.getCenter().lat, map.getCenter().lng, map._zoom ],
        "a": [ a.lat, a.lon ],
        "b": [ b.lat, b.lon ],
        "c": [ c.lat, c.lon ],
        "d": [ d.lat, d.lon ],
    };
    if (!autoStoring) {
        setInterval(function() {
            sessionStorage.setItem(stateName, JSON.stringify(init));
            sessionStorage.setItem(stateNamePosn, JSON.stringify(init.posn));
            sessionStorage.setItem("ellipsoidDisabled", isDisabled["ellipsoid"]);
            sessionStorage.setItem("sphereDisabled", isDisabled["sphere"]);
            
        }, 100);
        autoStoring = true;
    }
};

const map = L.map('map', { 
    zoomControl: false,
}).setView([init.posn[0], init.posn[1]], init.posn[2]);

L.tileLayer('https://tile.openstreetmap.org/{z}/{x}/{y}.png', {
    minZoom: 0,
    maxZoom: 17,
    preferCanvas: true,
    attribution: '&copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a>',
    noWrap: true,
}).addTo(map);

Module().then( Module => {

// document.getElementById('info-details').style.display = 'block';
if (opts.about) {
    document.getElementById('about').innerHTML = opts.about;
}

let shapesForDisabling = ['sphere', 'ellipsoid'];

for (let i = 0; i < shapesForDisabling.length; i++) {
    let shape = shapesForDisabling[i];

    if (isDisabled[shape]) {
        document.getElementById('toggle-'+shape).style.backgroundColor = 'inherit';
    }
    document.getElementById("toggle-"+shape).addEventListener("click", function() {
        if (isDisabled[shape]) {
            isDisabled[shape] = false;
            this.style.backgroundColor = '';    
        } else {
            isDisabled[shape] = true;
            this.style.backgroundColor = 'initial'
        }
        update();
    });
}


window.inverse = function(shape, lat1, lon1, lat2, lon2) {
    shape = shape == "ellipsoid" ? 1 : 0;
    var ptr = Module._malloc(24);
    Module.ccall("geokit_inverse", "null",
        ['number', 'number', 'number', 'number', 'number', 'number', 'number', 'number'],
        [shape, lat1, lon1, lat2, lon2, ptr, ptr+8, ptr+16]
    );
    var res = new Float64Array(Module.HEAPU8.buffer, ptr, 24);
    Module._free(ptr);
    return {
        s12: res[0],
        azi1: res[1],
        azi2: res[2],
    }
};

window.direct = function(shape, lat1, lon1, azi1, s12) {
    shape = shape == "ellipsoid" ? 1 : 0;
    var ptr = Module._malloc(24);
    Module.ccall("geokit_direct", "null",
        ['number', 'number', 'number', 'number', 'number', 'number', 'number', 'number'],
        [shape, lat1, lon1, azi1, s12, ptr, ptr+8, ptr+16]
    );
    var res = new Float64Array(Module.HEAPU8.buffer, ptr, 24);
    Module._free(ptr);
    return {
        lat2: res[0],
        lon2: res[1],
        azi2: res[2],
    }
}

window.track = function(shape, lat1, lon1, azi1, lat2, lon2) {
    shape = shape == "ellipsoid" ? 1 : 0;
    var ptr = Module._malloc(16);
    Module.ccall("geokit_track", "null",
        ['number', 'number', 'number', 'number', 'number', 'number', 'number', 'number'],
        [shape, lat1, lon1, azi1, lat2, lon2, ptr, ptr+8]
    );
    var res = new Float64Array(Module.HEAPU8.buffer, ptr, 16);
    Module._free(ptr);
    return {
        sat: res[0],
        sxt: res[1],
    };
}

window.intercept = function(shape, lat1, lon1, azi1, lat2, lon2, azi2) {
    shape = shape == "ellipsoid" ? 1 : 0;
    var ptr = Module._malloc(8*2);
    Module.ccall("geokit_intercept", "null",
        ['number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number'],
        [shape, lat1, lon1, azi1, lat2, lon2, azi2, ptr, ptr+8]
    );
    var res = new Float64Array(Module.HEAPU8.buffer, ptr, 16);
    Module._free(ptr);
    return {
        sat: res[0],
        sxt: res[1],
    }
};

window.intersect = function(shape, alat, alon, blat, blon, clat, clon, dlat, dlon) {
    shape = shape == "ellipsoid" ? 1 : 0;
    var ptr = Module._malloc(8*2);
    var isect = Module.ccall("geokit_intersect", 'number',
        ['number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number'],
        [shape, alat, alon, blat, blon, clat, clon, dlat, dlon, ptr, ptr+8]
    );
    var res = new Float64Array(Module.HEAPU8.buffer, ptr, 16);
    Module._free(ptr);
    return {
        isect: !!isect,
        lat: res[0],
        lon: res[1],
    }
};

window.nearest = function(shape, alat, alon, blat, blon, clat, clon, dlat, dlon) {
    shape = shape == "ellipsoid" ? 1 : 0;
    var ptr = Module._malloc(56);
    Module.ccall("geokit_nearest", 'null',
        ['number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 
         'number', 'number', 'number', 'number', 'number', 'number', 'number'],
        [shape, alat, alon, blat, blon, clat, clon, dlat, dlon, 
         ptr, ptr+8, ptr+16, ptr+24, ptr+32, ptr+40, ptr+48]
    );
    var res = new Float64Array(Module.HEAPU8.buffer, ptr, 56);
    Module._free(ptr);
    return {
        lat1: res[0],
        lon1: res[1],
        lat2: res[2],
        lon2: res[3],
        s12:  res[4],
        azi1: res[5],
        azi2: res[6],
    }
};

window.equibox = function(shape, alat, alon, blat, blon) {
    shape = shape == "ellipsoid" ? 1 : 0;
    var ptr = Module._malloc(32);
    Module.ccall("geokit_equibox", "null",
        ['number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number'],
        [shape, alat, alon, blat, blon, ptr, ptr+8, ptr+16, ptr+24]
    );
    var res = new Float64Array(Module.HEAPU8.buffer, ptr, 32);
    Module._free(ptr);
    return {
        minlat : res[0],
        minlon : res[1],
        maxlat : res[2],
        maxlon : res[3],
    }
};

window.equipath = function(shape, aminlat, aminlon, amaxlat, amaxlon, bminlat, bminlon, bmaxlat, bmaxlon) {
    shape = shape == "ellipsoid" ? 1 : 0;
    var ptr = Module._malloc(8*7);
    Module.ccall("geokit_equipath", "null",
        ['number', 'number', 'number', 'number', 'number', 'number', 'number', 'number',
         'number', 'number', 'number', 'number', 'number', 'number', 'number'],
        [shape, aminlat, aminlon, amaxlat, amaxlon, bminlat, bminlon, bmaxlat, bmaxlon, 
        ptr, ptr+8, ptr+16, ptr+24, ptr+32, ptr+40, ptr+48]
    );
    var res = new Float64Array(Module.HEAPU8.buffer, ptr, 56);
    Module._free(ptr);
    return {
        lat1 : res[0],
        lon1 : res[1],
        lat2 : res[2],
        lon2 : res[3],
        s12 : res[4],
        azi1 : res[5],
        azi2 : res[6],
    }
};

window.wrap360 = function(degs) {
    if (degs < 0 || degs >= 360) {
        degs = degs % 360.0;
        if (degs < 0) {
            degs += 360;
        } else if (degs >= 360) {
            degs -= 360;
        }
    }
    return degs;
}

// let reverseInverse = function(res) {
//     return {
//         s12: res.s12,
//         azi1: wrap360(res.azi2+180),
//         azi2: wrap360(res.azi1+180)
//     };
// }

function drawMapStraightLine(alat, alon, blat, blon, opts) {
    let opts2 = { weight: 3, color: "purple" };
    for (let name in opts) {
        opts2[name] = opts[name];
    }
    if (opts2.dashed && !opts2.dashArray) {
        opts2.dashArray = '1,10';
        delete opts2.dashed;
    }
    opts = opts2;
    return [L.polyline([[alat,alon],[blat,blon]], opts).addTo(map)];
    
}

function drawMapLine(shape, alat, alon, blat, blon, opts) {
    let opts2 = { weight: 3, color: "purple" };
    for (let name in opts) {
        opts2[name] = opts[name];
    }
    if (opts2.dashed && !opts2.dashArray) {
        opts2.dashArray = '1,10';
        delete opts2.dashed;
    }
    opts = opts2;
    let lines = [];
    let points = [];
    let plat, plon;
    function pushPoint(clat, clon) {
        if (!opts.noWrap && points.length && Math.abs(plon - clon) > 180) {
            if (plon - clon > 0) {
                // lift pen and move to other side of map
                points.push([clat, 180]);
                lines.push(L.polyline(points, opts).addTo(map));
                points = [[clat, -180]];
            } else {
                // lift pen and move to other side of map
                points.push([clat, -180]);
                lines.push(L.polyline(points, opts).addTo(map));
                points = [[clat, 180]];
            }
        }
        points.push([clat, clon]);
        plat = clat;
        plon = clon;
    }
    function pushPointSeries(alat, alon, azi1, s12) {
        let res;
        // One point every every 20km
        let npoints = Math.floor(s12/1000/20); 
        npoints = Math.min(Math.max(npoints, 10), 2000);
        for (let i = 0; i < npoints; i++) {
            res = direct(shape, alat, alon, azi1, (s12 / npoints) * (i+1));
            pushPoint(res.lat2, res.lon2);
        }
        return res;
    }
    let ab = inverse(shape, alat, alon, blat, blon);
    pushPoint(alat, alon);
    pushPointSeries(alat, alon, ab.azi1, ab.s12);
    let fline = L.polyline(points, opts).addTo(map)
    fline.command = ["drawMapLine", arguments];
    lines.push(fline);
    return lines;
}



window.toDist = function(x) {
    if (Math.abs(x) > 1000) {
        return (x/1000).toFixed(2) + " km";
    }
    if (Math.abs(x) > 1) {
        return x.toFixed(2) + " m";
    }
    return (x*1000).toFixed(2) + " cm";
}

function drawMapDot(lat1, lon1, opts) {
    let color = (opts||{}).color||"purple";
    let opts2 = { color: color, fillColor: color, weight: 2, radius: 6 };
    for (let name in opts) {
        opts2[name] = opts[name];
    }
    return L.circleMarker([lat1, lon1], opts2).addTo(map);
}

let getCapLatLng = function(marker, alt) {
    let ll = marker ? marker.getLatLng() : { lat: alt[0], lng: alt[1] };
    ll.lat = Math.max(-90, ll.lat);
    ll.lat = Math.min(90, ll.lat);
    ll.lng = Math.max(-179.9999, ll.lng);
    ll.lng = Math.min(179.9999, ll.lng);
    if (marker) { 
        marker.setLatLng(ll);
    }
    return ll;
}

let ctx = { shapes: {} };

function colorForShape(shape) {
    return shape == "ellipsoid" ? 
        // 'rgb(191, 70, 52)' // "red"
        "red"
        : 
        "blue";
        // 'rgb(85, 125, 193)' // "blue"
}

function drawAppend(key, shape, entry) {    
    if (Array.isArray(entry)) {
        for (let i = 0; i < entry.length; i++) {
            drawAppend(key, shape, entry[i]);
        }
        return;
    }
    if (!ctx.shapes) ctx.shapes = {};
    if (!ctx.shapes[shape]) ctx.shapes[shape] = {};
    if (!ctx.shapes[shape][key]) ctx.shapes[shape][key] = [];
    ctx.shapes[shape][key].push(entry);
}


window.drawClear = function(key, shape) {
    if (!ctx.shapes) return;
    if (!ctx.shapes[shape]) return;
    if (!ctx.shapes[shape][key]) return;
    for (let i = 0; i < ctx.shapes[shape][key].length; i++) {
        map.removeLayer(ctx.shapes[shape][key][i]);
    }
    delete ctx.shapes[shape][key];
}

window.drawLine = function(key, shape, lat1, lon1, lat2, lon2, opts) {
    opts = opts || {};
    opts.color = opts.color || colorForShape(shape);
    window.drawClear(key, shape);
    drawAppend(key, shape, drawMapLine(shape, lat1, lon1, lat2, lon2, opts));
};

window.drawDirectLine = function(key, shape, lat1, lon1, azi1, sat1, opts) {
	let final = direct(shape, lat1, lon1, azi1, sat1);
	const LIMIT = 5000000;
	let i = 0;
	let lat2 = lat1;
	let lon2 = lon1;
	let azi2 = azi1;
	let sat2 = sat1;
	while (sat2 > LIMIT) {
		let res = direct(shape, lat2, lon2, azi2, LIMIT);
		drawLine(key+":"+i, shape, lat2, lon2, res.lat2, res.lon2, opts);
		azi2 = res.azi2;
		lat2 = res.lat2;
		lon2 = res.lon2;
		sat2 -= LIMIT;	
		i++;
	}
	drawLine(key+":"+i, shape, lat2, lon2, final.lat2, final.lon2, opts);
}

window.drawStraightLine = function(key, shape, lat1, lon1, lat2, lon2, opts) {
    opts = opts || {};
    opts.color = opts.color || colorForShape(shape);
    window.drawClear(key, shape);
    drawAppend(key, shape, drawMapStraightLine(lat1, lon1, lat2, lon2, opts));
}



window.drawDot = function(key, shape, lat1, lon1, opts) {
    opts = opts || {};
    opts.color = opts.color || colorForShape(shape);
    window.drawClear(key, shape);
    drawAppend(key, shape, drawMapDot(lat1, lon1, opts));
};

window.mininfo = "";
let infolog = "";
window.info = {
    log: function() {
        let line = '';
        for (let i = 0; i < arguments.length; i++) {
            if (i > 0) line += ' ';
            line += arguments[i];
        }
        line += '\n';
        infolog += line;
    }
}

let markers = {};
function update() {
    var llA = getCapLatLng(markers['a'], init.a);
    var llB = getCapLatLng(markers['b'], init.b);
    var llC = getCapLatLng(markers['c'], init.c);
    var llD = getCapLatLng(markers['d'], init.d);
    let a = { lat: llA.lat, lon: llA.lng }; 
    let b = { lat: llB.lat, lon: llB.lng }; 
    let c = { lat: llC.lat, lon: llC.lng }; 
    let d = { lat: llD.lat, lon: llD.lng };
    if (opts.adjust) {
        let o = adjust(a, b, c, d);
        if (markers['a']) {
            markers['a'].setLatLng(new L.LatLng(o.a.lat, o.a.lon)); 
        }
        if (markers['b']) {
            markers['b'].setLatLng(new L.LatLng(o.b.lat, o.b.lon)); 
        }
        if (markers['c']) {
            markers['c'].setLatLng(new L.LatLng(o.c.lat, o.c.lon)); 
        }
        if (markers['d']) {
            markers['d'].setLatLng(new L.LatLng(o.d.lat, o.d.lon)); 
        }
    }

    saveState(a, b, c, d, map);
    
    let shapes = ["ellipsoid", "sphere"];
    for (let i = 0; i < shapes.length; i++) {
        let shape = shapes[i];
        if (ctx.shapes && ctx.shapes[shape]) {
            for (let key in ctx.shapes[shape]) {
                drawClear(key, shape);
            }
        }
        mininfo = "";
        infolog = "";
        if (!isDisabled[shape]) {
            opts.draw(shape, a, b, c, d, map.getZoomScale(1));
        }
        let el = document.getElementById('vals_'+shape)
        if (el) {
            el.innerText = infolog.trim();
        }
        el = document.getElementById('vals_'+shape+"_min")
        if (el) {
            el.innerText = mininfo.trim();
        }
    }
}

let toshow = ['a', 'b', 'c', 'd'];
for (let i = 0; i < toshow.length; i++) {
    let key = toshow[i];
    if (opts.show.includes(key)) {
        let name = (((opts.repl || {key:key})[key]) || key).toUpperCase();
        let marker = L.marker([init[key][0], init[key][1]], { draggable: 'true', 
            icon: L.divIcon({html: '<div><div>'+name+'</div></div>'})}).addTo(map);
        marker.on('drag', function() { update(); });
        markers[key] = marker;
    }
}

// let gcolor = '#98A';
// let gcolor = 'rgba(255, 240, 187, 0.1)';
let gcolor1 = 'rgba(64, 64, 64, 1)'; //'#AAA';
let gcolor2 = 'rgba(128, 128, 128, 1)'; //'#AAA';


// let gcolor1 = 'rgba(21, 19, 15)';
// let gcolor2 = 'rgba(21, 19, 15, 0.5)';

let gridBoundsOpts = { color: gcolor1, weight: 3.0 };
let gridMajorOpts = { color: gcolor1, weight: 1.5};
let gridMinorOpts = { color: gcolor1, weight: 0.25 };

if (1) {
    // grid lines
    const space = 30;
    drawMapStraightLine(90, 0, -90, 0, gridMajorOpts);
    drawMapStraightLine(0, -180, 0, 180, gridMajorOpts);
    // lat
    for (let lat = space; lat < 90; lat += space) {
        drawMapStraightLine(lat, -180, lat, 180, gridMinorOpts);    
    }
    for (let lat = -space; lat > -90; lat -= space) {
        drawMapStraightLine(lat, -180, lat, 180, gridMinorOpts);    
    }
    // lon
    for (let lon = space; lon < 180; lon += space) {
        drawMapStraightLine(90, lon, -90, lon, gridMinorOpts);
    }
    for (let lon = -space; lon > -180; lon -= space) {
        drawMapStraightLine(90, lon, -90, lon, gridMinorOpts);
    }
}

drawMapStraightLine(90, -180, 90, 180, gridBoundsOpts);
drawMapStraightLine(-90, -180, -90, 180, gridBoundsOpts);
drawMapStraightLine(90, 180, -90, 180, gridBoundsOpts);
drawMapStraightLine(90, -180, -90, -180, gridBoundsOpts);

map.on('zoom', update); 
map.on('move', update); 
update();


});
////////////////////////////////////////////////////////////////////////////////
}

