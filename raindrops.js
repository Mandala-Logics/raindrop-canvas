const cellCollisionMultiplier = 0.1;
const spreadingForceMultiplier = 0.5;
const polygonCollisionMultiplier = 0.01;
const polygonResolution = 8;
const g = 2;
const polygonSmoothLength = 10;
const smoothingIncriment = 0.5;
const smoothingRadius = 2;
const smoothLoopLimit = 1e3;
const maxCellCount = 400;
const cellRadius = 2;
const cellMax = 0.5;
const cellMin = 0.35;
const trailQ = 0.3;
const Tau = 2 * Math.PI;
const fuzziness = 0.4;
const drop_x_cohesion = cellRadius * 5;
const drop_y_cohesion = cellRadius * 10;
var RectCorner;
(function (RectCorner) {
    RectCorner[RectCorner["TopLeft"] = 0] = "TopLeft";
    RectCorner[RectCorner["TopRight"] = 1] = "TopRight";
    RectCorner[RectCorner["BottomLeft"] = 2] = "BottomLeft";
    RectCorner[RectCorner["BottomRight"] = 3] = "BottomRight";
})(RectCorner || (RectCorner = {}));
function DivRem(a, b) {
    var quotient = Math.trunc(a / b);
    var remainder = a - (quotient * b);
    return { quotient, remainder };
}
function GetBoundingRect(cells) {
    var x_max = -Infinity;
    var y_max = -Infinity;
    var x_min = Infinity;
    var y_min = Infinity;
    var top_left;
    var bottom_right;
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
export class RaindropCanvas {
    constructor(drawCallback) {
        this.Cells = [];
        this.Gravity = new PolarNumber(g, Math.PI / 2);
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
    }
    CellCount() { return this.Cells.length; }
    OnCanvasResize(HTMLCanvas) {
        var margin = 100;
        this.canvasRect = new Rect(new Vector(-margin, -margin), HTMLCanvas.width + margin, HTMLCanvas.height + margin);
    }
    ShedDrops() {
        if (this.Cells.length < maxCellCount) {
            this.Polygons.forEach(polygon => {
                if (Math.random() < 0.5) {
                    this.Cells.push(new DropCell(trailQ, polygon.Centroid));
                }
            });
        }
    }
    Stop() {
        this.Cells = [];
        this.keepLoopGoing = false;
    }
    Start() {
        requestAnimationFrame(this.Loop);
    }
    DefinePolygons() {
        this.Polygons = [];
        var clusters = ClusterCells(this.Cells, drop_x_cohesion, drop_y_cohesion);
        for (let i = 0; i < clusters.length; i++) {
            this.Polygons.push(new DropPolygon(clusters[i], polygonResolution));
        }
    }
    SpawnDrop(Q, x, y) {
        if (this.Cells.length < maxCellCount) {
            this.Cells.push(new DropCell(Q, new Vector(x, y)));
        }
    }
    ApplyPolygonCollision() {
        var circle;
        this.Polygons.forEach(polygon => {
            polygon.Cells.forEach(cell => {
                circle = cell.collision.GetShift(cell.force);
                var I = polygon.GetIntersection(circle);
                cell.ApplyForce(I, polygonCollisionMultiplier);
            });
        });
    }
    ApplyCellCollision() {
        var b;
        var i1;
        var i2;
        var c_a, c_b;
        var force;
        for (let a = 0; a < this.Cells.length; a++) {
            for (b = 0; b < this.Cells.length; b++) {
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
                    this.Cells[b].ApplyForce(force = new PolarNumber(i2.r - i1.r, i1.theta), cellCollisionMultiplier);
                    this.Cells[a].ApplyForce(force.GetInverse(), cellCollisionMultiplier);
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
        var rm = [];
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
        var r;
        var i;
        var add = [];
        this.Cells.forEach(cell => {
            if (cell.Q > cellMax) {
                r = DivRem(cell.Q, cellMax);
                cell.Q = Math.max(r.remainder, cellMin);
                for (i = 1; i <= r.quotient; i++) {
                    add.push(new DropCell(cellMax, cell.collision.center.Clone(0, 0)));
                }
            }
        });
        add.forEach(newCell => {
            this.Cells.push(newCell);
        });
    }
    ApplySpreadingForce() {
        var b, a;
        var z;
        this.Polygons.forEach(polygon => {
            for (a = 0; a < polygon.Cells.length; a++) {
                for (b = 0; b < polygon.Cells.length; b++) {
                    if (a === b) {
                        continue;
                    }
                    else {
                        z = polygon.Cells[a].collision.Intersection(polygon.Cells[b].collision);
                        if (z instanceof PolarNumber) {
                            polygon.Cells[b].ApplyForce(z, spreadingForceMultiplier);
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
class PolarNumber {
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
        var x = other.x - this.x;
        var y = other.y - this.y;
        return Math.atan2(y, x);
    }
    GetDistance(other) {
        var x = other.x - this.x;
        var y = other.y - this.y;
        return Math.sqrt(x * x + y * y);
    }
    GetMidpoint(other) {
        var dx = other.x - this.x;
        var dy = other.y - this.y;
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
    constructor(cells, pointsPerCircle = 16) {
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
    }
    TotalQ() {
        var ret = 0;
        this.Cells.forEach(cell => {
            ret += cell.Q;
        });
        return ret;
    }
    GetIntersection(circle) {
        var ret = new PolarNumber(0, 0);
        var x;
        this.Lines.forEach(line => {
            x = circle.GetLineIntersection(line);
            if (x instanceof PolarNumber) {
                ret.Add(x);
            }
        });
        return ret;
    }
    SmoothPolygon() {
        var smoothVectors = [];
        this.Lines.forEach(line => smoothVectors.push(...line.Break(polygonSmoothLength)));
        var circles = [];
        smoothVectors.forEach(vector => circles.push(new Circle(vector, smoothingRadius)));
        var changed = true;
        var z;
        var shifted;
        var c1, c2;
        var I;
        var locked = [];
        var n = 0;
        while (changed && n < smoothLoopLimit) {
            changed = false;
            n++;
            for (c1 = 0; c1 < circles.length; c1++) {
                if (locked.includes(circles[c1])) {
                    continue;
                }
                z = new PolarNumber(smoothingIncriment, (circles[c1].center.x > this.Centroid.x) ? Math.PI : 0);
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
        var hull = this.SmoothPolygon();
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
        var hull = this.SmoothPolygon();
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
    GetRandomInteralPoint() {
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
        context.moveTo(this.BottomRight.x, this.TopLeft.y);
        context.moveTo(this.BottomRight.x, this.BottomRight.y);
        context.moveTo(this.TopLeft.x, this.BottomRight.y);
        context.moveTo(this.TopLeft.x, this.TopLeft.y);
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
        context.arc(this.center.x, this.center.y, cellRadius, 0, Tau);
        context.stroke();
    }
    GetLineIntersection(line) {
        var nearestPoint;
        var d;
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
            let m = line.m();
            let c = line.c();
            let a = this.center.x;
            let b = this.center.y;
            let r = this.radius;
            var A = 1 + m * m;
            var B = 2 * m * (c - b) - 2 * a;
            var C = a * a + (c - b) * (c - b) - r * r;
            var D = B * B - 4 * A * C;
            if (D < 0) {
                return undefined;
            }
            else {
                var denominator = 2 * A;
                D = Math.sqrt(D);
                var x = (-B + D) / denominator;
                var p1 = new Vector(x, line.y(x));
                x = (-B - D) / denominator;
                nearestPoint = p1.GetMidpoint(new Vector(x, line.y(x)));
                d = nearestPoint.GetDistance(this.center);
            }
        }
        return new PolarNumber(d, Math.atan2(this.center.y - nearestPoint.y, this.center.x - nearestPoint.x));
    }
    CenterChordIntersects(m) {
        var denominator = Math.sqrt(1 + m * m);
        var a = this.radius / denominator;
        var b = (this.radius * m) / denominator;
        return [new Vector(this.center.x + a, this.center.y + b), new Vector(this.center.x - a, this.center.y - b)];
    }
    Intersection(other) {
        var dx = other.center.x - this.center.x;
        var dy = other.center.y - this.center.y;
        var r = Math.sqrt(dx * dx + dy * dy);
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
    constructor(Q, position) {
        this.force = new PolarNumber(0, 0);
        if (Q <= 0) {
            debugger;
        }
        this.Q = Q;
        this.collision = new Circle(position, cellRadius);
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
