//CONSTS
const Tau = 2 * Math.PI;
const fuzziness = 0.4;

enum RectCorner
{
    TopLeft = 0,
    TopRight = 1,
    BottomLeft = 2,
    BottomRight = 3
}

function DivRem(a: number, b: number) : { quotient: number; remainder: number }
{
    const quotient = Math.trunc(a / b);
    const remainder = a - (quotient * b);
    return { quotient, remainder };
}

function GetBoundingRect(cells : DropCell[]) : Rect
{
    let x_max : number = -Infinity;
    let y_max : number = -Infinity;
    let x_min : number = Infinity;
    let y_min : number = Infinity;
    let top_left : Vector;
    let bottom_right : Vector;

    cells.forEach(cell => {
        top_left = cell.collision.GetCorner(RectCorner.TopLeft);
        bottom_right = cell.collision.GetCorner(RectCorner.BottomRight);
        
        if (top_left.x < x_min) { x_min = top_left.x; }

        if (top_left.y < y_min) { y_min = top_left.y }

        if (bottom_right.x > x_max) { x_max = bottom_right.x; }

        if (bottom_right.y > y_max) { y_max = bottom_right.y; }
    }); 

    return new Rect(new Vector(x_min, y_min), x_max - x_min, y_max - y_min);
}

function ClusterCells(cells: DropCell[], max_x: number, max_y: number): DropCell[][] {
    const clusters: DropCell[][] = [];
    const visited = new Set<DropCell>();

    for (const start of cells) {
        if (visited.has(start)) continue;

        const cluster: DropCell[] = [];
        const stack: DropCell[] = [start];

        while (stack.length > 0) {
            const current = stack.pop()!;
            if (visited.has(current)) continue;

            visited.add(current);
            cluster.push(current);

            for (const other of cells) {
                if (visited.has(other)) continue;
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
  // MULTIPLIERS
  public cellCollisionMultiplier: number;
  public spreadingForceMultiplier: number;
  public polygonCollisionMultiplier: number;
  public polygonResolution: number;
  public gravity: number;

  // SMOOTHING
  public polygonSmoothLength: number;
  public smoothingIncriment: number;
  public smoothingRadius: number;
  public smoothLoopLimit: number;

  // MAFS
  public maxCellCount: number;
  public cellRadius: number;
  public cellMax: number;
  public cellMin: number;
  public trailQ: number;
  public drop_x_cohesion: number;
  public drop_y_cohesion: number;
  public gravity_direction : number;

  //MARGIN
  public canvasMargin : number;

  constructor(overrides: Partial<RaindropSettings> = {}) {
    // MULTIPLIERS
    this.cellCollisionMultiplier = overrides.cellCollisionMultiplier ?? 0.1;
    this.spreadingForceMultiplier = overrides.spreadingForceMultiplier ?? 0.5;
    this.polygonCollisionMultiplier = overrides.polygonCollisionMultiplier ?? 0.01;
    this.polygonResolution = overrides.polygonResolution ?? 8;
    this.gravity = overrides.gravity ?? 2;

    // SMOOTHING
    this.polygonSmoothLength = overrides.polygonSmoothLength ?? 10;
    this.smoothingIncriment = overrides.smoothingIncriment ?? 0.5;
    this.smoothingRadius = overrides.smoothingRadius ?? 2;
    this.smoothLoopLimit = overrides.smoothLoopLimit ?? 1e3;

    // MAFS
    this.maxCellCount = overrides.maxCellCount ?? 400;
    this.cellRadius = overrides.cellRadius ?? 2;
    this.cellMax = overrides.cellMax ?? 0.5;
    this.cellMin = overrides.cellMin ?? 0.35;
    this.trailQ = overrides.trailQ ?? 0.3;
    this.drop_x_cohesion = overrides.drop_x_cohesion ?? (this.cellRadius * 5);
    this.drop_y_cohesion = overrides.drop_y_cohesion ?? (this.cellRadius * 10);
    this.gravity_direction = overrides.gravity_direction ?? Math.PI / 2;

    //MARGIN
    this.canvasMargin = overrides.canvasMargin ?? 50;
  }
}

export class RaindropCanvas
{
    private Cells : Array<DropCell> = []
    private Gravity : PolarNumber;
    private Polygons : Array<DropPolygon> = []
    private canvasRect : Rect = new Rect(new Vector(0, 0), 0, 0);
    private keepLoopGoing : boolean = true;
    private callback : (polygons : DropPolygon[]) => void;
    private settings : RaindropSettings;

    constructor(drawCallback : (polygons : DropPolygon[]) => void)
    {
        this.callback = drawCallback;
        this.settings = new RaindropSettings();

        this.Gravity = new PolarNumber(this.settings.gravity, this.settings.gravity_direction);
    }

    public ApplySettings(settings: Partial<RaindropSettings>) {
        Object.assign(this.settings, settings);
        
        this.Gravity = new PolarNumber(this.settings.gravity, this.settings.gravity_direction);
    }

    public GetSettings(): RaindropSettings {
        return this.settings;
    }

    CellCount() { return this.Cells.length; }

    OnCanvasResize(HTMLCanvas : HTMLCanvasElement)
    {
        const margin = this.settings.canvasMargin;

        this.canvasRect = new Rect(new Vector(-margin, -margin), HTMLCanvas.width + margin, HTMLCanvas.height + margin);
    }

    ShedDrops()
    {
        if (this.Cells.length < this.settings.maxCellCount)
        {
            this.Polygons.forEach(polygon => {
                if (Math.random() < 0.5)
                {
                    this.Cells.push(new DropCell(this.settings.trailQ, polygon.Centroid, this.settings.cellRadius));
                }
            });
        }
    }

    Stop() : void
    {
        this.Cells = [];
        this.Polygons = [];
        this.keepLoopGoing = false;
    }

    Start()
    {
        requestAnimationFrame(this.Loop);
    }

    Loop = (currentTime: DOMHighResTimeStamp): void => {

        if (this.Cells.length > 0)
        {
            this.RemoveCells();

            this.SplitCells();

            this.ApplySpreadingForce();

            this.ApplyCellCollision();

            this.ApplyPolygonCollision();

            this.MoveCells();

            this.DefinePolygons();
        }

        this.callback(this.Polygons);

        if (this.keepLoopGoing) { requestAnimationFrame(this.Loop); } // keep the loop going
    }

    DefinePolygons() : void
    {
        this.Polygons = [];

        const clusters = ClusterCells(this.Cells, this.settings.drop_x_cohesion, this.settings.drop_y_cohesion);

        for (let i = 0; i < clusters.length; i++)
        {
            this.Polygons.push(new DropPolygon(this.settings, clusters[i], this.settings.polygonResolution));
        }
    }

    SpawnDrop(Q : number, x : number, y : number): void
    {
        if (this.Cells.length < this.settings.maxCellCount)
        {
            this.Cells.push(new DropCell(Q, new Vector(x, y), this.settings.cellRadius));
        }        
    }

    ApplyPolygonCollision() : void
    {
        this.Polygons.forEach(polygon => {
            polygon.Cells.forEach(cell => {

                let circle = cell.collision.GetShift(cell.force);

                let I = polygon.GetIntersection(circle);

                cell.ApplyForce(I, this.settings.polygonCollisionMultiplier);

            });
        });
    }

    ApplyCellCollision() : void
    {
        let i1 : PolarNumber | undefined;
        let i2 : PolarNumber | undefined;
        let c_a : Circle, c_b : Circle;
        let force : PolarNumber;

        for (let a = 0; a < this.Cells.length; a++)
        {
            for (let b = 0; b < this.Cells.length; b++)
            {
                if (a === b) { continue; }

                c_a = this.Cells[a].collision;
                c_b = this.Cells[b].collision;

                if (c_a.IsStacked(c_b)) { continue; }
                
                i1 = c_a.Intersection(c_b)

                if (!i1) { continue; }

                i2 = c_a.GetShift(this.Cells[a].force).Intersection(c_b.GetShift(this.Cells[b].force));

                if (!i2) { continue; }

                if (i2.r > i1.r) //a moving into b
                {
                    this.Cells[b].ApplyForce(force = new PolarNumber(i2.r - i1.r, i1.theta), this.settings.cellCollisionMultiplier);
                    this.Cells[a].ApplyForce(force.GetInverse(), this.settings.cellCollisionMultiplier);
                }
            }
        }
    }

    MoveCells()
    {
        this.Cells.forEach(cell => {
            cell.Shift(cell.force);
            cell.force = this.Gravity.Clone(cell.Q);
        });
    }

    RemoveCells()
    {
        const rm : Array<DropCell> = []

        this.Cells.forEach(cell => {
            if (!this.canvasRect.IsContained(cell.collision.center))
            {
                rm.push(cell);
            }
        });

        if (rm.length > 0)
        {
            this.Cells = this.Cells.filter((cell) => { return !rm.includes(cell); });
        }
    }

    SplitCells() : void
    {
        let r : { quotient: number; remainder: number } ;
        let i : number;
        const add : Array<DropCell> = []

        this.Cells.forEach(cell => {

            if (cell.Q > this.settings.cellMax)
            {
                r = DivRem(cell.Q, this.settings.cellMax);

                cell.Q = Math.max(r.remainder, this.settings.cellMin);

                for (i = 1; i <= r.quotient; i++)
                {
                    add.push(new DropCell(this.settings.cellMax, cell.collision.center.Clone(0, 0), this.settings.cellRadius));
                }
            }
            
        });

        add.forEach(newCell => {
            this.Cells.push(newCell);
        });
    }

    ApplySpreadingForce() : void
    {
        let b : number, a : number;
        let z : PolarNumber | undefined;

        this.Polygons.forEach(polygon => {

            for (a = 0; a < polygon.Cells.length; a++)
            {
                for (b = 0; b < polygon.Cells.length; b++)
                {
                    if (a === b) { continue; }
                    else
                    {
                        z = polygon.Cells[a].collision.Intersection(polygon.Cells[b].collision);
    
                        if (z instanceof PolarNumber) //overlapping or stacked
                        {
                            polygon.Cells[b].ApplyForce(z, this.settings.spreadingForceMultiplier);
                        }
                    }                
                }
            }

        });
    }

    ApplyForce(force : PolarNumber) : void
    {
        this.Cells.forEach(cell => {
            cell.ApplyForce(force)
        });
    }
}

export class PolarNumber
{
    public r : number;
    public theta : number;

    constructor(r: number, theta : number)
    {
        this.r = r;
        this.theta = theta;
    }

    Add(other: PolarNumber, factor : number = 1) : void
    {        
        const x1 = this.r * Math.cos(this.theta);
        const y1 = this.r * Math.sin(this.theta);
    
        const x2 = other.r * factor * Math.cos(other.theta);
        const y2 = other.r * factor * Math.sin(other.theta);
        
        const x = x1 + x2;
        const y = y1 + y2;
        
        this.r = Math.sqrt(x * x + y * y);
        this.theta = Math.atan2(y, x);

        if (this.theta < 0) { this.theta += Tau; }
    }   
    
    Clone(factor : number = 1) : PolarNumber
    {
        return new PolarNumber(this.r * factor, this.theta);
    }

    GetInverse() : PolarNumber
    {
        return new PolarNumber(this.r, this.theta + Math.PI)
    }
}

class Vector
{
    public x : number
    public y : number

    constructor(x: number, y : number)
    {
        this.x = x;
        this.y = y;
    }

    Shift(shift: PolarNumber): void {
        const dx = shift.r * Math.cos(shift.theta);
        const dy = shift.r * Math.sin(shift.theta);
    
        this.x += dx;
        this.y += dy;
    }

    GetShift(shift: PolarNumber) : Vector {
        const dx = shift.r * Math.cos(shift.theta);
        const dy = shift.r * Math.sin(shift.theta);
    
        return new Vector(this.x + dx, this.y + dy)
    }

    GetDirection(other : Vector) : number
    {
        const x = other.x - this.x;
        const y = other.y - this.y;

        return Math.atan2(y, x)
    }

    GetDistance(other : Vector) : number
    {
        const x = other.x - this.x;
        const y = other.y - this.y;

        return Math.sqrt(x * x + y * y);
    }

    GetMidpoint(other : Vector) : Vector
    {
        const dx = other.x - this.x;
        const dy = other.y - this.y;

        return new Vector(this.x + dx / 2, this.y + dy / 2);
    }

    Clone(dx : number, dy : number) : Vector
    {
        return new Vector(this.x + dx, this.y + dy);
    }
}

class Line
{
    public v1 : Vector;
    public v2 : Vector;

    constructor(v1 : Vector, v2 : Vector)
    {
        this.v1 = v1;
        this.v2 = v2;
    }

    length() : number { return this.v1.GetDistance(this.v2); }

    Break(length: number): Vector[] {
        const totalLength = this.length();
        const r = DivRem(totalLength, length);
    
        if (r.quotient > 0) {
            const segments = r.quotient;
            const step = 1 / (segments + 1); // number of *gaps* between points
            const ret: Vector[] = [this.v1];
    
            for (let i = 1; i <= segments; i++) {
                const t = step * i;
                const x = this.v1.x + (this.v2.x - this.v1.x) * t;
                const y = this.v1.y + (this.v2.y - this.v1.y) * t;
                ret.push(new Vector(x, y));
            }
    
            ret.push(this.v2);
            return ret;
        } else {
            return [this.v1, this.v2];
        }
    }
    

    m() : number
    {
        if (this.IsHorizontal()) { return Infinity; }
        else { return (this.v2.x - this.v1.x)/(this.v2.y - this.v1.y) }
    }

    c() : number | undefined
    {
        if (this.IsVertical()) { return undefined; }
        else { return this.v1.y - this.m() * this.v1.x; }
    }

    y(x : number) : number | undefined
    {
        if (this.IsVertical()) { return undefined; }
        else if (this.IsHorizontal()) { return this.v1.y; }
        else { return this.m() * x + (this.c() as number); }
    }

    x(y : number) : number | undefined
    {
        if (this.IsVertical()) { return this.v1.x; }
        else if (this.IsHorizontal()) { return undefined }
        else { return (y - (this.c() as number)) / this.m(); }        
    }

    x_max() : number { return Math.max(this.v1.x, this.v2.x); }

    x_min() : number { return Math.min(this.v1.x, this.v2.x); }

    IsVertical() : boolean { return this.v2.x - this.v1.x < fuzziness; }

    IsHorizontal() : boolean { return this.v2.y - this.v1.y < fuzziness; }
}

export class DropPolygon {
    public Hull: Vector[] = [];
    public Cells : DropCell[];
    public Lines : Line[] = [];
    public Centroid : Vector;
    private settings : RaindropSettings;

    constructor(settings : RaindropSettings, cells: DropCell[], pointsPerCircle: number = 16) {
        const allPoints: Vector[] = [];

        // Sample edge points from each circle
        for (const cell of cells) {
            allPoints.push(...DropPolygon.sampleCircle(cell, pointsPerCircle));
        }

        // Compute convex hull
        this.Hull = DropPolygon.convexHull(allPoints);
        this.Cells = cells;

        this.GetLines();
        this.Centroid = GetBoundingRect(this.Cells).GetCenter();

        this.settings = settings;
    }

    TotalQ() : number
    {
        let ret = 0;

        this.Cells.forEach(cell => {
            ret += cell.Q;
        });

        return ret;
    }

    public GetIntersection(circle : Circle) : PolarNumber
    {
        const ret = new PolarNumber(0, 0);
        let x : PolarNumber | undefined;

        this.Lines.forEach(line => {
            
            x = circle.GetLineIntersection(line);

            if (x instanceof PolarNumber)
            {
                ret.Add(x);
            }

        });

        return ret;
    }

    SmoothPolygon() : Vector[]
    {
        const smoothVectors : Vector[] = [];

        this.Lines.forEach(line => smoothVectors.push(...line.Break(this.settings.polygonSmoothLength)));

        const circles : Circle[] = [];

        smoothVectors.forEach(vector => circles.push(new Circle(vector, this.settings.smoothingRadius)));

        let changed : boolean = true;
        let z : PolarNumber;
        let shifted : Circle;
        let c1 : number, c2 : number;
        let I : PolarNumber | undefined;
        const locked : Circle[] = []
        let n = 0;

        while (changed && n < this.settings.smoothLoopLimit)
        {
            changed = false;
            n++;

            for (c1 = 0; c1 < circles.length; c1++)
            {
                if (locked.includes(circles[c1])) { continue; }

                z = new PolarNumber(this.settings.smoothingIncriment, (circles[c1].center.x > this.Centroid.x) ? Math.PI : 0);

                shifted = circles[c1].GetShift(z);

                for (c2 = 0; c2 < circles.length; c2++)
                {
                    if (c1 === c2) { continue; }

                    I = shifted.Intersection(circles[c2]);

                    if (I instanceof PolarNumber)
                    {
                        locked.push(circles[c1]);
                        break;
                    }
                }

                if (I instanceof PolarNumber) { continue; }

                for (c2 = 0; c2 < this.Cells.length; c2++)
                {
                    I = shifted.Intersection(this.Cells[c2].collision);

                    if (I instanceof PolarNumber)
                    {
                        locked.push(circles[c1]);
                        break;
                    }
                }

                if (!I)
                {
                    changed = true;
                    circles[c1].center = smoothVectors[c1] = shifted.center;
                }
            }
        }

        return smoothVectors;   
    }

    GetLines() : void
    {
        for (let i = 0; i < this.Hull.length - 1; i++) {
            this.Lines.push(new Line(this.Hull[i], this.Hull[i + 1]));
        }

        this.Lines.push(new Line(this.Hull[this.Hull.length - 1], this.Hull[0]));
    }

    private static sampleCircle(cell: DropCell, steps: number): Vector[] {
        const points: Vector[] = [];
        for (let i = 0; i < steps; i++) {
            const angle = (2 * Math.PI * i) / steps;
            points.push(new Vector(
                cell.collision.center.x + cell.collision.radius * Math.cos(angle),
                cell.collision.center.y + cell.collision.radius * Math.sin(angle)
            ));
        }
        return points;
    }

    private static cross(o: Vector, a: Vector, b: Vector): number {
        return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
    }

    private static convexHull(points: Vector[]): Vector[] {
        const sorted = [...points].sort((a, b) =>
            a.x === b.x ? a.y - b.y : a.x - b.x
        );

        const lower: Vector[] = [];
        for (const p of sorted) {
            while (lower.length >= 2 && DropPolygon.cross(
                lower[lower.length - 2],
                lower[lower.length - 1],
                p
            ) <= 0) {
                lower.pop();
            }
            lower.push(p);
        }

        const upper: Vector[] = [];
        for (let i = sorted.length - 1; i >= 0; i--) {
            const p = sorted[i];
            while (upper.length >= 2 && DropPolygon.cross(
                upper[upper.length - 2],
                upper[upper.length - 1],
                p
            ) <= 0) {
                upper.pop();
            }
            upper.push(p);
        }

        // Remove the last point of each half (they're duplicates)
        upper.pop();
        lower.pop();

        return lower.concat(upper);
    }

    public Draw(context: CanvasRenderingContext2D): void {
        if (this.Hull.length === 0) return;

        context.beginPath();
        context.moveTo(this.Hull[0].x, this.Hull[0].y);
        for (let i = 1; i < this.Hull.length; i++) {
            context.lineTo(this.Hull[i].x, this.Hull[i].y);
        }
        context.closePath();
        context.stroke();
    }

    public DrawSmooth(context: CanvasRenderingContext2D): void {
        if (this.Hull.length === 0) return;

        const hull = this.SmoothPolygon();

        context.beginPath();
        context.moveTo(hull[0].x, hull[0].y);
        for (let i = 1; i < hull.length; i++) {
            context.lineTo(hull[i].x, hull[i].y);
        }
        context.closePath();
        context.stroke();
    }

    public Fill(context: CanvasRenderingContext2D): void {
        if (this.Hull.length === 0) return;

        context.beginPath();
        context.moveTo(this.Hull[0].x, this.Hull[0].y);
        for (let i = 1; i < this.Hull.length; i++) {
            context.lineTo(this.Hull[i].x, this.Hull[i].y);
        }
        context.closePath();
        context.fill();
    }

    public FillSmooth(context: CanvasRenderingContext2D): void {
        if (this.Hull.length === 0) return;

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

class Rect
{
    public TopLeft : Vector;
    public BottomRight : Vector;

    constructor(topLeft : Vector, width : number, height : number)
    {
        this.TopLeft = topLeft.Clone(0,0);
        this.BottomRight = new Vector(topLeft.x + width, topLeft.y + height)
    }

    GetRandomInternalPoint() : Vector
    {
        return new Vector(this.Width() * Math.random() + this.TopLeft.x, this.Height() * Math.random() + this.TopLeft.y);
    }

    IsContained(vector : Vector) : boolean
    {
        return vector.x >= this.TopLeft.x
            && vector.x <= this.BottomRight.x
            && vector.y >= this.TopLeft.y
            && vector.y <= this.BottomRight.y;
    }

    Width() : number { return this.BottomRight.x - this.TopLeft.x; }

    Height() : number { return this.BottomRight.y - this.TopLeft.y; }

    GetCenter() : Vector { return this.TopLeft.GetMidpoint(this.BottomRight); }

    GetCorner(corner : RectCorner) : Vector
    {
        switch (corner)
        {
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

    Draw(context : CanvasRenderingContext2D) : void {
        context.beginPath();
        context.moveTo(this.TopLeft.x, this.TopLeft.y);
        context.lineTo(this.BottomRight.x, this.TopLeft.y);
        context.lineTo(this.BottomRight.x, this.BottomRight.y);
        context.lineTo(this.TopLeft.x, this.BottomRight.y);
        context.closePath();
        context.stroke();
    }
}

class Circle
{
    public center : Vector;
    public radius : number;

    constructor(center : Vector, radius : number)
    {
        this.center = center;
        this.radius = radius;
    }

    IsStacked(other : Circle) : boolean
    {
        return this.center.x - other.center.x < fuzziness
            && this.center.y - other.center.y < fuzziness;
    }

    GetShift(shift : PolarNumber) : Circle
    {
        return new Circle(this.center.GetShift(shift), this.radius);
    }

    GetCorner(corner : RectCorner)
    {
        switch (corner)
        {
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

    Width() : number { return this.radius * 2; }

    Height() : number { return this.radius * 2; }

    Draw(context : CanvasRenderingContext2D)
    {
        context.beginPath();
        context.arc(this.center.x, this.center.y, this.radius, 0, Tau);
        context.stroke();
    }

    GetLineIntersection(line : Line) : PolarNumber | undefined
    {
        let nearestPoint : Vector;
        let d : number;

        if (line.IsHorizontal())
        {
            if (line.v1.y <= this.center.y + this.radius && line.v1.y >= this.center.y - this.radius)
            {
                nearestPoint = new Vector(this.center.x, line.v1.y);
                d = Math.abs(line.v1.y - this.center.y);
            }
            else { return undefined; }
        }
        else if (line.IsVertical())
        {
            if (line.v1.x <= this.center.x + this.radius && line.v1.x >= this.center.x - this.radius)
            {
                nearestPoint = new Vector(line.v1.x, this.center.y);
                d = Math.abs(line.v1.x - this.center.x);
            }
            else { return undefined; }
        }
        else
        {
            const m = line.m();
            const c = line.c() as number;
            const a = this.center.x;
            const b = this.center.y
            const r = this.radius;
            
            const A = 1 + m * m;
            const B = 2 * m * (c - b) - 2 * a
            const C = a * a + (c - b) * (c - b) - r * r;

            let D = B * B - 4 * A * C;

            if (D < 0) { return undefined; }
            else
            {
                const denominator = 2 * A;
                D = Math.sqrt(D);

                let x = (-B + D) / denominator;

                const p1 = new Vector(x, line.y(x) as number);

                x = (-B - D) / denominator;

                nearestPoint = p1.GetMidpoint(new Vector(x, line.y(x) as number));

                d = nearestPoint.GetDistance(this.center);
            }
        }
        
        return new PolarNumber(d, Math.atan2(this.center.y - nearestPoint.y, this.center.x - nearestPoint.x));
    }

    CenterChordIntersects(m : number) : [Vector, Vector]
    {
        const denominator = Math.sqrt(1 + m * m);

        const a = this.radius / denominator;
        const b = (this.radius * m) / denominator;

        return [new Vector(this.center.x + a, this.center.y + b), new Vector(this.center.x - a, this.center.y - b)]
    }

    Intersection(other : Circle) : PolarNumber | undefined
    {
        const dx = other.center.x - this.center.x;
        const dy = other.center.y - this.center.y;

        const r = Math.sqrt(dx * dx + dy * dy);

        if (this.radius + other.radius - r < fuzziness) { return undefined; } //don't touch
        else if (r < fuzziness) { return new PolarNumber(this.radius + other.radius, Math.random() * Tau); } //stacked
        else { return new PolarNumber(this.radius + other.radius - r, Math.atan2(dy, dx)); } //overlapping
    }
}

class DropCell
{
    public Q : number;
    public collision : Circle;
    public force : PolarNumber = new PolarNumber(0, 0);

    constructor(Q : number, position : Vector, radius : number)
    {
        if (Q <= 0) { debugger; }

        this.Q = Q;
        this.collision = new Circle(position, radius);
    }

    ApplyForce(force : PolarNumber, scale : number = 1) : void
    {
        this.force.Add(force, scale * this.Q);
    }

    Shift(x : PolarNumber)
    {
        this.collision.center.Shift(x);
    }

    Draw(context : CanvasRenderingContext2D)
    {
        this.collision.Draw(context);
    }
}

