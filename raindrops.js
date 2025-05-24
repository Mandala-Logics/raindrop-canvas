/*!
 * Raindrop Canvas Effect
 * 
 * © Mandala Logics – https://operalimina.co.uk
 * GitHub: https://github.com/Mandala-Logics/raindrop-canvas
 * Buy Me a Coffee: https://www.buymeacoffee.com/operalimina
 * 
 * This code is free to use for non-commercial purposes.
 * If you use it, please credit operalimina.co.uk.
 */

const Tau = 2 * Math.PI;
const fuzziness = 0.4;
var RectCorner;
(function (RectCorner) {
    RectCorner[RectCorner["TopLeft"] = 0] = "TopLeft";
    RectCorner[RectCorner["TopRight"] = 1] = "TopRight";
    RectCorner[RectCorner["BottomLeft"] = 2] = "BottomLeft";
    RectCorner[RectCorner["BottomRight"] = 3] = "BottomRight";
})(RectCorner || (RectCorner = {}));
function DivRem(a, b) {
    const quotient = Math.trunc(a / b);
    const remainder = a - (quotient * b);
    return { quotient, remainder };
}
function GetBoundingRect(cells) {
    let x_max = -Infinity;
    let y_max = -Infinity;
    let x_min = Infinity;
    let y_min = Infinity;
    let top_left;
    let bottom_right;
    cells.forEach(cell => {
        top_left = cell.collision.GetCorner(RectCorner.TopLeft);
        bottom_right = cell.collision.GetCorner(RectCorner.BottomRight);
        if (top_left.x < x_min) {
            x_min = top_left.x;
        }
        if (top_left.y < y_min) {
            y_min = top_left.y;
        }
        if (bottom_right.x > x_max) {
            x_max = bottom_right.x;
        }
        if (bottom_right.y > y_max) {
            y_max = bottom_right.y;
        }
    });
    return new Rect(new Vector(x_min, y_min), x_max - x_min, y_max - y_min);
}
function ClusterCells(cells, max_x, max_y) {
    const clusters = [];
    const visited = new Set();
    for (const start of cells) {
        if (visited.has(start))
            continue;
        const cluster = [];
        const stack = [start];
        while (stack.length > 0) {
            const current = stack.pop();
            if (visited.has(current))
                continue;
            visited.add(current);
            cluster.push(current);
            for (const other of cells) {
                if (visited.has(other))
                    continue;
                const dx = Math.abs(current.collision.center.x - other.collision.center.x);
                const dy = Math.abs(current.collision.center.y - other.collision.center.y);
                if (dx <= max_x && dy <= max_y) {
                    stack.push(other);
                }
            }
        }
        clusters.push(cluster);
    }
    return clusters;
}
export class RaindropSettings {
    constructor(overrides = {}) {
        var _a, _b, _c, _d, _e, _f, _g, _h, _j, _k, _l, _m, _o, _p, _q, _r, _s, _t;
        this.cellCollisionMultiplier = (_a = overrides.cellCollisionMultiplier) !== null && _a !== void 0 ? _a : 0.1;
        this.spreadingForceMultiplier = (_b = overrides.spreadingForceMultiplier) !== null && _b !== void 0 ? _b : 0.5;
        this.polygonCollisionMultiplier = (_c = overrides.polygonCollisionMultiplier) !== null && _c !== void 0 ? _c : 0.01;
        this.polygonResolution = (_d = overrides.polygonResolution) !== null && _d !== void 0 ? _d : 8;
        this.gravity = (_e = overrides.gravity) !== null && _e !== void 0 ? _e : 2;
        this.polygonSmoothLength = (_f = overrides.polygonSmoothLength) !== null && _f !== void 0 ? _f : 10;
        this.smoothingIncriment = (_g = overrides.smoothingIncriment) !== null && _g !== void 0 ? _g : 0.5;
        this.smoothingRadius = (_h = overrides.smoothingRadius) !== null && _h !== void 0 ? _h : 2;
        this.smoothLoopLimit = (_j = overrides.smoothLoopLimit) !== null && _j !== void 0 ? _j : 1e3;
        this.maxCellCount = (_k = overrides.maxCellCount) !== null && _k !== void 0 ? _k : 400;
        this.cellRadius = (_l = overrides.cellRadius) !== null && _l !== void 0 ? _l : 2;
        this.cellMax = (_m = overrides.cellMax) !== null && _m !== void 0 ? _m : 0.5;
        this.cellMin = (_o = overrides.cellMin) !== null && _o !== void 0 ? _o : 0.35;
        this.trailQ = (_p = overrides.trailQ) !== null && _p !== void 0 ? _p : 0.3;
        this.drop_x_cohesion = (_q = overrides.drop_x_cohesion) !== null && _q !== void 0 ? _q : (this.cellRadius * 5);
        this.drop_y_cohesion = (_r = overrides.drop_y_cohesion) !== null && _r !== void 0 ? _r : (this.cellRadius * 10);
        this.gravity_direction = (_s = overrides.gravity_direction) !== null && _s !== void 0 ? _s : Math.PI / 2;
        this.canvasMargin = (_t = overrides.canvasMargin) !== null && _t !== void 0 ? _t : 50;
    }
}
export class RaindropCanvas {
    constructor(drawCallback) {
        this.Cells = [];
        this.Polygons = [];
        this.canvasRect = new Rect(new Vector(0, 0), 0, 0);
        this.keepLoopGoing = true;
        this.Loop = (currentTime) => {
            if (this.Cells.length > 0) {
                this.RemoveCells();
                this.SplitCells();
                this.ApplySpreadingForce();
                this.ApplyCellCollision();
                this.ApplyPolygonCollision();
                this.MoveCells();
                this.DefinePolygons();
            }
            this.callback(this.Polygons);
            if (this.keepLoopGoing) {
                requestAnimationFrame(this.Loop);
            }
        };
        this.callback = drawCallback;
        this.settings = new RaindropSettings();
        this.Gravity = new PolarNumber(this.settings.gravity, this.settings.gravity_direction);
    }
    ApplySettings(settings) {
        Object.assign(this.settings, settings);
        this.Gravity = new PolarNumber(this.settings.gravity, this.settings.gravity_direction);
    }
    GetSettings() {
        return this.settings;
    }
    CellCount() { return this.Cells.length; }
    OnCanvasResize(HTMLCanvas) {
        const margin = this.settings.canvasMargin;
        this.canvasRect = new Rect(new Vector(-margin, -margin), HTMLCanvas.width + margin, HTMLCanvas.height + margin);
    }
    ShedDrops() {
        if (this.Cells.length < this.settings.maxCellCount) {
            this.Polygons.forEach(polygon => {
                if (Math.random() < 0.5) {
                    this.Cells.push(new DropCell(this.settings.trailQ, polygon.Centroid, this.settings.cellRadius));
                }
            });
        }
    }
    Stop() {
        this.Cells = [];
        this.Polygons = [];
        this.keepLoopGoing = false;
    }
    Start() {
        requestAnimationFrame(this.Loop);
    }
    DefinePolygons() {
        this.Polygons = [];
        const clusters = ClusterCells(this.Cells, this.settings.drop_x_cohesion, this.settings.drop_y_cohesion);
        for (let i = 0; i < clusters.length; i++) {
            this.Polygons.push(new DropPolygon(this.settings, clusters[i], this.settings.polygonResolution));
        }
    }
    SpawnDrop(Q, x, y) {
        if (this.Cells.length < this.settings.maxCellCount) {
            this.Cells.push(new DropCell(Q, new Vector(x, y), this.settings.cellRadius));
        }
    }
    ApplyPolygonCollision() {
        this.Polygons.forEach(polygon => {
            polygon.Cells.forEach(cell => {
                let circle = cell.collision.GetShift(cell.force);
                let I = polygon.GetIntersection(circle);
                cell.ApplyForce(I, this.settings.polygonCollisionMultiplier);
            });
        });
    }
    ApplyCellCollision() {
        let i1;
        let i2;
        let c_a, c_b;
        let force;
        for (let a = 0; a < this.Cells.length; a++) {
            for (let b = 0; b < this.Cells.length; b++) {
                if (a === b) {
                    continue;
                }
                c_a = this.Cells[a].collision;
                c_b = this.Cells[b].collision;
                if (c_a.IsStacked(c_b)) {
                    continue;
                }
                i1 = c_a.Intersection(c_b);
                if (!i1) {
                    continue;
                }
                i2 = c_a.GetShift(this.Cells[a].force).Intersection(c_b.GetShift(this.Cells[b].force));
                if (!i2) {
                    continue;
                }
                if (i2.r > i1.r) {
                    this.Cells[b].ApplyForce(force = new PolarNumber(i2.r - i1.r, i1.theta), this.settings.cellCollisionMultiplier);
                    this.Cells[a].ApplyForce(force.GetInverse(), this.settings.cellCollisionMultiplier);
                }
            }
        }
    }
    MoveCells() {
        this.Cells.forEach(cell => {
            cell.Shift(cell.force);
            cell.force = this.Gravity.Clone(cell.Q);
        });
    }
    RemoveCells() {
        const rm = [];
        this.Cells.forEach(cell => {
            if (!this.canvasRect.IsContained(cell.collision.center)) {
                rm.push(cell);
            }
        });
        if (rm.length > 0) {
            this.Cells = this.Cells.filter((cell) => { return !rm.includes(cell); });
        }
    }
    SplitCells() {
        let r;
        let i;
        const add = [];
        this.Cells.forEach(cell => {
            if (cell.Q > this.settings.cellMax) {
                r = DivRem(cell.Q, this.settings.cellMax);
                cell.Q = Math.max(r.remainder, this.settings.cellMin);
                for (i = 1; i <= r.quotient; i++) {
                    add.push(new DropCell(this.settings.cellMax, cell.collision.center.Clone(0, 0), this.settings.cellRadius));
                }
            }
        });
        add.forEach(newCell => {
            this.Cells.push(newCell);
        });
    }
    ApplySpreadingForce() {
        let b, a;
        let z;
        this.Polygons.forEach(polygon => {
            for (a = 0; a < polygon.Cells.length; a++) {
                for (b = 0; b < polygon.Cells.length; b++) {
                    if (a === b) {
                        continue;
                    }
                    else {
                        z = polygon.Cells[a].collision.Intersection(polygon.Cells[b].collision);
                        if (z instanceof PolarNumber) {
                            polygon.Cells[b].ApplyForce(z, this.settings.spreadingForceMultiplier);
                        }
                    }
                }
            }
        });
    }
    ApplyForce(force) {
        this.Cells.forEach(cell => {
            cell.ApplyForce(force);
        });
    }
}
export class PolarNumber {
    constructor(r, theta) {
        this.r = r;
        this.theta = theta;
    }
    Add(other, factor = 1) {
        const x1 = this.r * Math.cos(this.theta);
        const y1 = this.r * Math.sin(this.theta);
        const x2 = other.r * factor * Math.cos(other.theta);
        const y2 = other.r * factor * Math.sin(other.theta);
        const x = x1 + x2;
        const y = y1 + y2;
        this.r = Math.sqrt(x * x + y * y);
        this.theta = Math.atan2(y, x);
        if (this.theta < 0) {
            this.theta += Tau;
        }
    }
    Clone(factor = 1) {
        return new PolarNumber(this.r * factor, this.theta);
    }
    GetInverse() {
        return new PolarNumber(this.r, this.theta + Math.PI);
    }
}
class Vector {
    constructor(x, y) {
        this.x = x;
        this.y = y;
    }
    Shift(shift) {
        const dx = shift.r * Math.cos(shift.theta);
        const dy = shift.r * Math.sin(shift.theta);
        this.x += dx;
        this.y += dy;
    }
    GetShift(shift) {
        const dx = shift.r * Math.cos(shift.theta);
        const dy = shift.r * Math.sin(shift.theta);
        return new Vector(this.x + dx, this.y + dy);
    }
    GetDirection(other) {
        const x = other.x - this.x;
        const y = other.y - this.y;
        return Math.atan2(y, x);
    }
    GetDistance(other) {
        const x = other.x - this.x;
        const y = other.y - this.y;
        return Math.sqrt(x * x + y * y);
    }
    GetMidpoint(other) {
        const dx = other.x - this.x;
        const dy = other.y - this.y;
        return new Vector(this.x + dx / 2, this.y + dy / 2);
    }
    Clone(dx, dy) {
        return new Vector(this.x + dx, this.y + dy);
    }
}
class Line {
    constructor(v1, v2) {
        this.v1 = v1;
        this.v2 = v2;
    }
    length() { return this.v1.GetDistance(this.v2); }
    Break(length) {
        const totalLength = this.length();
        const r = DivRem(totalLength, length);
        if (r.quotient > 0) {
            const segments = r.quotient;
            const step = 1 / (segments + 1);
            const ret = [this.v1];
            for (let i = 1; i <= segments; i++) {
                const t = step * i;
                const x = this.v1.x + (this.v2.x - this.v1.x) * t;
                const y = this.v1.y + (this.v2.y - this.v1.y) * t;
                ret.push(new Vector(x, y));
            }
            ret.push(this.v2);
            return ret;
        }
        else {
            return [this.v1, this.v2];
        }
    }
    m() {
        if (this.IsHorizontal()) {
            return Infinity;
        }
        else {
            return (this.v2.x - this.v1.x) / (this.v2.y - this.v1.y);
        }
    }
    c() {
        if (this.IsVertical()) {
            return undefined;
        }
        else {
            return this.v1.y - this.m() * this.v1.x;
        }
    }
    y(x) {
        if (this.IsVertical()) {
            return undefined;
        }
        else if (this.IsHorizontal()) {
            return this.v1.y;
        }
        else {
            return this.m() * x + this.c();
        }
    }
    x(y) {
        if (this.IsVertical()) {
            return this.v1.x;
        }
        else if (this.IsHorizontal()) {
            return undefined;
        }
        else {
            return (y - this.c()) / this.m();
        }
    }
    x_max() { return Math.max(this.v1.x, this.v2.x); }
    x_min() { return Math.min(this.v1.x, this.v2.x); }
    IsVertical() { return this.v2.x - this.v1.x < fuzziness; }
    IsHorizontal() { return this.v2.y - this.v1.y < fuzziness; }
}
export class DropPolygon {
    constructor(settings, cells, pointsPerCircle = 16) {
        this.Hull = [];
        this.Lines = [];
        const allPoints = [];
        for (const cell of cells) {
            allPoints.push(...DropPolygon.sampleCircle(cell, pointsPerCircle));
        }
        this.Hull = DropPolygon.convexHull(allPoints);
        this.Cells = cells;
        this.GetLines();
        this.Centroid = GetBoundingRect(this.Cells).GetCenter();
        this.settings = settings;
    }
    TotalQ() {
        let ret = 0;
        this.Cells.forEach(cell => {
            ret += cell.Q;
        });
        return ret;
    }
    GetIntersection(circle) {
        const ret = new PolarNumber(0, 0);
        let x;
        this.Lines.forEach(line => {
            x = circle.GetLineIntersection(line);
            if (x instanceof PolarNumber) {
                ret.Add(x);
            }
        });
        return ret;
    }
    SmoothPolygon() {
        const smoothVectors = [];
        this.Lines.forEach(line => smoothVectors.push(...line.Break(this.settings.polygonSmoothLength)));
        const circles = [];
        smoothVectors.forEach(vector => circles.push(new Circle(vector, this.settings.smoothingRadius)));
        let changed = true;
        let z;
        let shifted;
        let c1, c2;
        let I;
        const locked = [];
        let n = 0;
        while (changed && n < this.settings.smoothLoopLimit) {
            changed = false;
            n++;
            for (c1 = 0; c1 < circles.length; c1++) {
                if (locked.includes(circles[c1])) {
                    continue;
                }
                z = new PolarNumber(this.settings.smoothingIncriment, (circles[c1].center.x > this.Centroid.x) ? Math.PI : 0);
                shifted = circles[c1].GetShift(z);
                for (c2 = 0; c2 < circles.length; c2++) {
                    if (c1 === c2) {
                        continue;
                    }
                    I = shifted.Intersection(circles[c2]);
                    if (I instanceof PolarNumber) {
                        locked.push(circles[c1]);
                        break;
                    }
                }
                if (I instanceof PolarNumber) {
                    continue;
                }
                for (c2 = 0; c2 < this.Cells.length; c2++) {
                    I = shifted.Intersection(this.Cells[c2].collision);
                    if (I instanceof PolarNumber) {
                        locked.push(circles[c1]);
                        break;
                    }
                }
                if (!I) {
                    changed = true;
                    circles[c1].center = smoothVectors[c1] = shifted.center;
                }
            }
        }
        return smoothVectors;
    }
    GetLines() {
        for (let i = 0; i < this.Hull.length - 1; i++) {
            this.Lines.push(new Line(this.Hull[i], this.Hull[i + 1]));
        }
        this.Lines.push(new Line(this.Hull[this.Hull.length - 1], this.Hull[0]));
    }
    static sampleCircle(cell, steps) {
        const points = [];
        for (let i = 0; i < steps; i++) {
            const angle = (2 * Math.PI * i) / steps;
            points.push(new Vector(cell.collision.center.x + cell.collision.radius * Math.cos(angle), cell.collision.center.y + cell.collision.radius * Math.sin(angle)));
        }
        return points;
    }
    static cross(o, a, b) {
        return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
    }
    static convexHull(points) {
        const sorted = [...points].sort((a, b) => a.x === b.x ? a.y - b.y : a.x - b.x);
        const lower = [];
        for (const p of sorted) {
            while (lower.length >= 2 && DropPolygon.cross(lower[lower.length - 2], lower[lower.length - 1], p) <= 0) {
                lower.pop();
            }
            lower.push(p);
        }
        const upper = [];
        for (let i = sorted.length - 1; i >= 0; i--) {
            const p = sorted[i];
            while (upper.length >= 2 && DropPolygon.cross(upper[upper.length - 2], upper[upper.length - 1], p) <= 0) {
                upper.pop();
            }
            upper.push(p);
        }
        upper.pop();
        lower.pop();
        return lower.concat(upper);
    }
    Draw(context) {
        if (this.Hull.length === 0)
            return;
        context.beginPath();
        context.moveTo(this.Hull[0].x, this.Hull[0].y);
        for (let i = 1; i < this.Hull.length; i++) {
            context.lineTo(this.Hull[i].x, this.Hull[i].y);
        }
        context.closePath();
        context.stroke();
    }
    DrawSmooth(context) {
        if (this.Hull.length === 0)
            return;
        const hull = this.SmoothPolygon();
        context.beginPath();
        context.moveTo(hull[0].x, hull[0].y);
        for (let i = 1; i < hull.length; i++) {
            context.lineTo(hull[i].x, hull[i].y);
        }
        context.closePath();
        context.stroke();
    }
    Fill(context) {
        if (this.Hull.length === 0)
            return;
        context.beginPath();
        context.moveTo(this.Hull[0].x, this.Hull[0].y);
        for (let i = 1; i < this.Hull.length; i++) {
            context.lineTo(this.Hull[i].x, this.Hull[i].y);
        }
        context.closePath();
        context.fill();
    }
    FillSmooth(context) {
        if (this.Hull.length === 0)
            return;
        const hull = this.SmoothPolygon();
        context.beginPath();
        context.moveTo(hull[0].x, hull[0].y);
        for (let i = 1; i < hull.length; i++) {
            context.lineTo(hull[i].x, hull[i].y);
        }
        context.closePath();
        context.fill();
    }
}
class Rect {
    constructor(topLeft, width, height) {
        this.TopLeft = topLeft.Clone(0, 0);
        this.BottomRight = new Vector(topLeft.x + width, topLeft.y + height);
    }
    GetRandomInternalPoint() {
        return new Vector(this.Width() * Math.random() + this.TopLeft.x, this.Height() * Math.random() + this.TopLeft.y);
    }
    IsContained(vector) {
        return vector.x >= this.TopLeft.x
            && vector.x <= this.BottomRight.x
            && vector.y >= this.TopLeft.y
            && vector.y <= this.BottomRight.y;
    }
    Width() { return this.BottomRight.x - this.TopLeft.x; }
    Height() { return this.BottomRight.y - this.TopLeft.y; }
    GetCenter() { return this.TopLeft.GetMidpoint(this.BottomRight); }
    GetCorner(corner) {
        switch (corner) {
            case RectCorner.TopLeft:
                return this.TopLeft;
            case RectCorner.TopRight:
                return new Vector(this.BottomRight.x, this.TopLeft.y);
            case RectCorner.BottomLeft:
                return new Vector(this.TopLeft.x, this.BottomRight.y);
            case RectCorner.BottomRight:
                return this.BottomRight;
        }
    }
    Draw(context) {
        context.beginPath();
        context.moveTo(this.TopLeft.x, this.TopLeft.y);
        context.lineTo(this.BottomRight.x, this.TopLeft.y);
        context.lineTo(this.BottomRight.x, this.BottomRight.y);
        context.lineTo(this.TopLeft.x, this.BottomRight.y);
        context.closePath();
        context.stroke();
    }
}
class Circle {
    constructor(center, radius) {
        this.center = center;
        this.radius = radius;
    }
    IsStacked(other) {
        return this.center.x - other.center.x < fuzziness
            && this.center.y - other.center.y < fuzziness;
    }
    GetShift(shift) {
        return new Circle(this.center.GetShift(shift), this.radius);
    }
    GetCorner(corner) {
        switch (corner) {
            case RectCorner.TopLeft:
                return new Vector(this.center.x - this.radius, this.center.y - this.radius);
            case RectCorner.TopRight:
                return new Vector(this.center.x + this.radius, this.center.y - this.radius);
            case RectCorner.BottomLeft:
                return new Vector(this.center.x - this.radius, this.center.y + this.radius);
            case RectCorner.BottomRight:
                return new Vector(this.center.x + this.radius, this.center.y + this.radius);
        }
    }
    Width() { return this.radius * 2; }
    Height() { return this.radius * 2; }
    Draw(context) {
        context.beginPath();
        context.arc(this.center.x, this.center.y, this.radius, 0, Tau);
        context.stroke();
    }
    GetLineIntersection(line) {
        let nearestPoint;
        let d;
        if (line.IsHorizontal()) {
            if (line.v1.y <= this.center.y + this.radius && line.v1.y >= this.center.y - this.radius) {
                nearestPoint = new Vector(this.center.x, line.v1.y);
                d = Math.abs(line.v1.y - this.center.y);
            }
            else {
                return undefined;
            }
        }
        else if (line.IsVertical()) {
            if (line.v1.x <= this.center.x + this.radius && line.v1.x >= this.center.x - this.radius) {
                nearestPoint = new Vector(line.v1.x, this.center.y);
                d = Math.abs(line.v1.x - this.center.x);
            }
            else {
                return undefined;
            }
        }
        else {
            const m = line.m();
            const c = line.c();
            const a = this.center.x;
            const b = this.center.y;
            const r = this.radius;
            const A = 1 + m * m;
            const B = 2 * m * (c - b) - 2 * a;
            const C = a * a + (c - b) * (c - b) - r * r;
            let D = B * B - 4 * A * C;
            if (D < 0) {
                return undefined;
            }
            else {
                const denominator = 2 * A;
                D = Math.sqrt(D);
                let x = (-B + D) / denominator;
                const p1 = new Vector(x, line.y(x));
                x = (-B - D) / denominator;
                nearestPoint = p1.GetMidpoint(new Vector(x, line.y(x)));
                d = nearestPoint.GetDistance(this.center);
            }
        }
        return new PolarNumber(d, Math.atan2(this.center.y - nearestPoint.y, this.center.x - nearestPoint.x));
    }
    CenterChordIntersects(m) {
        const denominator = Math.sqrt(1 + m * m);
        const a = this.radius / denominator;
        const b = (this.radius * m) / denominator;
        return [new Vector(this.center.x + a, this.center.y + b), new Vector(this.center.x - a, this.center.y - b)];
    }
    Intersection(other) {
        const dx = other.center.x - this.center.x;
        const dy = other.center.y - this.center.y;
        const r = Math.sqrt(dx * dx + dy * dy);
        if (this.radius + other.radius - r < fuzziness) {
            return undefined;
        }
        else if (r < fuzziness) {
            return new PolarNumber(this.radius + other.radius, Math.random() * Tau);
        }
        else {
            return new PolarNumber(this.radius + other.radius - r, Math.atan2(dy, dx));
        }
    }
}
class DropCell {
    constructor(Q, position, radius) {
        this.force = new PolarNumber(0, 0);
        if (Q <= 0) {
            debugger;
        }
        this.Q = Q;
        this.collision = new Circle(position, radius);
    }
    ApplyForce(force, scale = 1) {
        this.force.Add(force, scale * this.Q);
    }
    Shift(x) {
        this.collision.center.Shift(x);
    }
    Draw(context) {
        this.collision.Draw(context);
    }
}
