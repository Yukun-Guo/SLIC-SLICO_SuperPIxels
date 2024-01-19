function rgbToLab(rin:any, gin:any, bin:any, sz:any, lvec: any[], avec: any[], bvec: any[]) {
    let epsilon = 0.008856;
    let kappa = 903.3;
    
    let Xr = 0.950456;
    let Yr = 1.0;
    let Zr = 1.088754;
    
    for (let i = 0; i < sz; i++) {
        let sR = rin[i] / 255.0;
        let sG = gin[i] / 255.0;
        let sB = bin[i] / 255.0;

        let R = (sR <= 0.04045) ? sR / 12.92 : Math.pow((sR + 0.055) / 1.055, 2.4);
        let G = (sG <= 0.04045) ? sG / 12.92 : Math.pow((sG + 0.055) / 1.055, 2.4);
        let B = (sB <= 0.04045) ? sB / 12.92 : Math.pow((sB + 0.055) / 1.055, 2.4);

        let X = R * 0.4124564 + G * 0.3575761 + B * 0.1804375;
        let Y = R * 0.2126729 + G * 0.7151522 + B * 0.0721750;
        let Z = R * 0.0193339 + G * 0.1191920 + B * 0.9503041;

        let xr = X / Xr;
        let yr = Y / Yr;
        let zr = Z / Zr;

        let epsilonCheck = (val: number) => (val > epsilon) ? Math.pow(val, 1.0 / 3.0) : (kappa * val + 16.0) / 116.0;

        let fx = epsilonCheck(xr);
        let fy = epsilonCheck(yr);
        let fz = epsilonCheck(zr);

        let lval = 116.0 * fy - 16.0;
        let aval = 500.0 * (fx - fy);
        let bval = 200.0 * (fy - fz);

        lvec[i] = lval;
        avec[i] = aval;
        bvec[i] = bval;
    }
}


function getLABXYSeeds(STEP: number, width: number, height: number, seedIndices: any[], numseeds: number | undefined) {
    let hexgrid = false;
    let xstrips = Math.ceil(width / STEP);
    let ystrips = Math.ceil(height / STEP);

    let xerr = width - STEP * xstrips;
    if (xerr < 0) {
        xstrips--;
        xerr = width - STEP * xstrips;
    }

    let yerr = height - STEP * ystrips;
    if (yerr < 0) {
        ystrips--;
        yerr = height - STEP * ystrips;
    }

    let xerrperstrip = xerr / xstrips;
    let yerrperstrip = yerr / ystrips;

    let xoff = STEP / 2;
    let yoff = STEP / 2;

    let n = 0;
    for (let y = 0; y < ystrips; y++) {
        let ye = y * yerrperstrip;
        for (let x = 0; x < xstrips; x++) {
            let xe = x * xerrperstrip;
            let seedx = x * STEP + xoff + xe;
            
            if (hexgrid) {
                seedx = x * STEP + (xoff << (y & 0x1)) + xe;
                if (seedx >= width) seedx = width - 1;
            }

            let seedy = y * STEP + yoff + ye;
            let i = seedy * width + seedx;
            seedIndices[n] = i;
            n++;
        }
    }
    
    numseeds = n;
}

function performSuperpixelSLICO(lvec: any[], avec: any[], bvec: any[], kseedsl: any[], kseedsa: any[], kseedsb: any[], kseedsx: any[], kseedsy: any[], width: number, height: number, numseeds: number, klabels: any[], STEP: number) {
    let dx4 = [-1, 0, 1, 0];
    let dy4 = [0, -1, 0, 1];
    let sz = width * height;
    let numk = numseeds;
    let offset = STEP;

    let clustersize = new Array(numk).fill(0);
    let inv = new Array(numk).fill(0);
    let sigmal = new Array(numk).fill(0);
    let sigmaa = new Array(numk).fill(0);
    let sigmab = new Array(numk).fill(0);
    let sigmax = new Array(numk).fill(0);
    let sigmay = new Array(numk).fill(0);
    let distvec = new Array(sz).fill(Number.MAX_VALUE);
    let distlab = new Array(sz).fill(Number.MAX_VALUE);
    let maxlab = new Array(numk).fill(100.0 * 100.0);

    let invxywt = 1.0 / (STEP * STEP);

    for (let itr = 0; itr < 10; itr++) {
        for (let i = 0; i < sz; i++) {
            distvec[i] = Number.MAX_VALUE;
        }

        for (let n = 0; n < numk; n++) {
            let x1 = kseedsx[n] - offset;
            let y1 = kseedsy[n] - offset;
            let x2 = kseedsx[n] + offset;
            let y2 = kseedsy[n] + offset;

            for (let y = y1; y < y2; y++) {
                for (let x = x1; x < x2; x++) {
                    let i = y * width + x;
                    let l = lvec[i];
                    let a = avec[i];
                    let b = bvec[i];

                    distlab[i] = (l - kseedsl[n]) * (l - kseedsl[n]) +
                        (a - kseedsa[n]) * (a - kseedsa[n]) +
                        (b - kseedsb[n]) * (b - kseedsb[n]);

                    let distxy = (x - kseedsx[n]) * (x - kseedsx[n]) +
                        (y - kseedsy[n]) * (y - kseedsy[n]);

                    let dist = distlab[i] / maxlab[n] + distxy * invxywt;

                    if (dist < distvec[i]) {
                        distvec[i] = dist;
                        klabels[i] = n;
                    }
                }
            }
        }

        if (itr === 0) {
            for (let n = 0; n < numk; n++) {
                maxlab[n] = 1.0;
            }
        }

        for (let i = 0; i < sz; i++) {
            if (maxlab[klabels[i]] < distlab[i]) {
                maxlab[klabels[i]] = distlab[i];
            }
        }

        for (let k = 0; k < numk; k++) {
            sigmal[k] = 0;
            sigmaa[k] = 0;
            sigmab[k] = 0;
            sigmax[k] = 0;
            sigmay[k] = 0;
            clustersize[k] = 0;
        }

        let ind = 0;
        for (let r = 0; r < height; r++) {
            for (let c = 0; c < width; c++) {
                if (klabels[ind] >= 0) {
                    sigmal[klabels[ind]] += lvec[ind];
                    sigmaa[klabels[ind]] += avec[ind];
                    sigmab[klabels[ind]] += bvec[ind];
                    sigmax[klabels[ind]] += c;
                    sigmay[klabels[ind]] += r;
                    clustersize[klabels[ind]] += 1.0;
                }
                ind++;
            }
        }

        for (let k = 0; k < numk; k++) {
            if (clustersize[k] <= 0) {
                clustersize[k] = 1;
            }
            inv[k] = 1.0 / clustersize[k];
        }

        for (let k = 0; k < numk; k++) {
            kseedsl[k] = sigmal[k] * inv[k];
            kseedsa[k] = sigmaa[k] * inv[k];
            kseedsb[k] = sigmab[k] * inv[k];
            kseedsx[k] = sigmax[k] * inv[k];
            kseedsy[k] = sigmay[k] * inv[k];
        }
    }
}


function enforceSuperpixelConnectivity(labels:any, width:any, height:any, numSuperpixels:any, nlabels:any, finalNumberOfLabels:any) {
    let dx4 = [-1, 0, 1, 0];
    let dy4 = [0, -1, 0, 1];
    let sz = width * height;
    let SUPSZ = Math.ceil(sz / numSuperpixels);
    let xvec = new Array(SUPSZ * 10);
    let yvec = new Array(SUPSZ * 10);

    for (let i = 0; i < sz; i++) {
        nlabels[i] = -1;
    }

    let oindex = 0;
    let adjlabel = 0;
    let label = 0;

    for (let j = 0; j < height; j++) {
        for (let k = 0; k < width; k++) {
            if (nlabels[oindex] < 0) {
                nlabels[oindex] = label;

                let x1 = k;
                let y1 = j;
                let x2 = k + 1;
                let y2 = j + 1;

                for (let n = 0; n < 4; n++) {
                    let x = x1 + dx4[n];
                    let y = y1 + dy4[n];
                    if ((x >= 0 && x < width) && (y >= 0 && y < height)) {
                        let nindex = y * width + x;
                        if (nlabels[nindex] >= 0) adjlabel = nlabels[nindex];
                    }
                }

                let count = 1;
                for (let c = 0; c < count; c++) {
                    for (let n = 0; n < 4; n++) {
                        let x = xvec[c] + dx4[n];
                        let y = yvec[c] + dy4[n];

                        if ((x >= 0 && x < width) && (y >= 0 && y < height)) {
                            let nindex = y * width + x;

                            if (nlabels[nindex] < 0 && labels[oindex] === labels[nindex]) {
                                xvec[count] = x;
                                yvec[count] = y;
                                nlabels[nindex] = label;
                                count++;
                            }
                        }
                    }
                }

                if (count <= SUPSZ >> 2) {
                    for (let c = 0; c < count; c++) {
                        let ind = yvec[c] * width + xvec[c];
                        nlabels[ind] = adjlabel;
                    }
                    label--;
                }
                label++;
            }
            oindex++;
        }
    }
    finalNumberOfLabels = label;

    return {
        nlabels: nlabels,
        finalNumberOfLabels: finalNumberOfLabels
    };
}


function SLICO(imgData: { data: any; width: any; height: any; }, options: { regionSize: any; minRegionSize?: number; }) {
    let height = imgData.height;
    let width = imgData.width;
    let numRegionsX = width / options.regionSize;//@ts-ignore
    let numRegionsY = height / options.regionSize;//@ts-ignore
    let numSuperpixels = numRegionsX * numRegionsY;//@ts-ignore
    let sz = width * height;

    // reformat image data from array of pixels to 2d array of rgb values
    let img = new Array(height);
    for (let y = 0; y < height; y++) {
        img[y] = new Array(width);
        for (let x = 0; x < width; x++) {
            img[y][x] = imgData.data.slice((y * width + x) * 4, (y * width + x) * 4 + 3);
        }
    }

    let rin = new Array(sz);
    let gin = new Array(sz);
    let bin = new Array(sz);

    let lvec = new Array(sz);
    let avec = new Array(sz);
    let bvec = new Array(sz);

    // Convert RGB to LAB
    for (let x = 0, ii = 0; x < width; x++) {
        for (let y = 0; y < height; y++) {
            let i = y * width + x;
            rin[i] = img[y][x][0];
            gin[i] = img[y][x][1];
            bin[i] = img[y][x][2];
        }
    }

    rgbToLab(rin, gin, bin, sz, lvec, avec, bvec);

    // Get LAB seeds
    let step = Math.sqrt(sz / numSuperpixels) + 0.5;
    let seedIndices = new Array(sz);
    let numseeds;
    getLABXYSeeds(step, width, height, seedIndices, numseeds);

    let kseedsx = new Array(numseeds);
    let kseedsy = new Array(numseeds);
    let kseedsl = new Array(numseeds);
    let kseedsa = new Array(numseeds);
    let kseedsb = new Array(numseeds);
    //@ts-ignore
    for (let k = 0; k < numseeds; k++) { 
        kseedsx[k] = seedIndices[k] % width;
        kseedsy[k] = Math.floor(seedIndices[k] / width);
        kseedsl[k] = lvec[seedIndices[k]];
        kseedsa[k] = avec[seedIndices[k]];
        kseedsb[k] = bvec[seedIndices[k]];
    }

    // Perform superpixel SLICO
    let klabels = new Array(sz);//@ts-ignore
    performSuperpixelSLICO(lvec, avec, bvec, kseedsl, kseedsa, kseedsb, kseedsx, kseedsy, width, height, numseeds, klabels, step);

    // Enforce connectivity
    let nlabels = new Array(sz);
    let finalNumberOfLabels;
    enforceSuperpixelConnectivity(klabels, width, height, numSuperpixels, nlabels, finalNumberOfLabels);

    return nlabels;
}


// Public API.
export default function (img: Uint8Array, width: number, height: number, pixelSize: number = 40, minPixelSize?: number) {
    let options = { regionSize: pixelSize, minRegionSize: minPixelSize || (pixelSize * pixelSize) / 4 };
  
    let imageData = {
      data: img,
      width: width,
      height: height
    };
    var segmentation = SLICO(imageData, options);
    return segmentation;
  }