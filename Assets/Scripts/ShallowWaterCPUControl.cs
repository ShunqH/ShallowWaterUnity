using System;
using System.Drawing;
using System.Net.Http;
using Unity.VisualScripting;
using UnityEngine;
using UnityEngine.Rendering;

//[ExecuteAlways]
public class ShallowWater_CPUSim : MonoBehaviour
{
    [Header("Shallow Water Simulation Compute Shader")]
    public ComputeShader SWSimShader; 
    [Header("Render Material")]
    public Material WaterMaterial;
    public Material SeaBedMaterial;

    [Header("General")]
    public float g = 1.0f;
    public float H0 = 1.0f;
    [Header("Simulation")]
    public float CFL = 0.5f;
    public float maxTimePerFrame = 0.0166667f;
    public int maxStepPerFrame = 100; 
    public float injectWaveFreq = 0.05f;
    public float injectWaveAmp = 0.2f; 
    [Header("Grid")]
    public float xmin = 0f;
    public float xmax = 10f;
    public float zmin = -5f;
    public float zmax = 5f;
    public int NX = 1024; 
    public int NZ = 1024;
    [Header("Sea Bed")]
    public float bmin = 0.0f;
    public float bmax = 2.0f;
    public float bAmp = 0.2f;
    public float bZWaveNumber = 3f; 

    private float hmin = 1e-6f;
    private float dx, dz, drmin; 
    private int NGhost = 1;
    private int lenx, lenz, lenArr;
    private float cmax;
    private float[] x, z; 
    private float[] h, velx, velz, b, eta; 
    private Vector3[] fHalf; 
    private Vector3[] gHalf; 
    private float[] parbx, parbz;

    private float tnow = 0f;
    private float dt = 0.016667f;
    private float dtdx, dtdz;

    private int THREADS_X = 16;
    private int THREADS_Z = 16;
    private int dispatchX, dispatchZ, dispatch_BDX_X, dispatch_BDX_Z, dispatch_BDZ_X, dispatch_BDZ_Z;
    private int kernelUpdateBDX, kernelUpdateBDZ, kernelUpdateF, kernelUpdateG, kernelUpdateU;

    private ComputeBuffer hBuffer, velxBuffer, velzBuffer, bBuffer, etaBuffer;
    private ComputeBuffer fHalfBuffer, gHalfBuffer;
    private ComputeBuffer parbxBuffer, parbzBuffer;

    private void Awake()
    {
        lenx = NX + NGhost * 2;
        lenz = NZ + NGhost * 2;
        lenArr = lenx * lenz;
        CreateBuffer();
        float[] init = new float[lenArr];
        for (int i = 0; i < lenArr; ++i)
            init[i] = H0;

        etaBuffer.SetData(init);

        WaterMaterial.SetBuffer("eta", etaBuffer);
    }

    // Start is called once before the first execution of Update after the MonoBehaviour is created
    void Start()
    {
        Setup();
        BuiltMesh();
        FindKernel();
        WaterMaterial.SetBuffer("eta", etaBuffer);
        //CreateBuffer();
        bindBuffer();
        SetBufferData();
        SetCBuffer();       // all constants.
        
    }

    // Update is called once per frame
    void Update()
    {
        float tThisFrame = Mathf.Min(Time.deltaTime, maxTimePerFrame);
        tThisFrame = Mathf.Max(tThisFrame, 1e-5f);
        dt = CFL * drmin / cmax;
        int step = Mathf.CeilToInt(tThisFrame / dt);
        dt = tThisFrame / step;
        step = Mathf.Min(maxStepPerFrame, step); 

        dtdx = dt / dx;
        dtdz = dt / dz;
        SWSimShader.SetFloat("dt", dt);
        SWSimShader.SetFloat("dtdx", dtdx);
        SWSimShader.SetFloat("dtdz", dtdz);
        for (int nt = 0; nt < step; nt++) { 
            tnow += dt;
            SWSimShader.SetFloat("tnow", tnow);

            SWSimShader.Dispatch(kernelUpdateBDX, dispatch_BDX_X, dispatch_BDX_Z, 1);
            SWSimShader.Dispatch(kernelUpdateBDZ, dispatch_BDZ_X, dispatch_BDZ_Z, 1);
            SWSimShader.Dispatch(kernelUpdateF, dispatchX, dispatchZ, 1);
            SWSimShader.Dispatch(kernelUpdateG, dispatchX, dispatchZ, 1);
            SWSimShader.Dispatch(kernelUpdateU, dispatchX, dispatchZ, 1);
            WaterMaterial.SetBuffer("eta", etaBuffer);
        }
        //Debug.Log($"Step per frame: {step}");
    }
    
    void Setup()
    {
        lenx = NX + NGhost * 2;
        lenz = NZ + NGhost * 2;
        lenArr = lenx * lenz;
        cmax = Mathf.Sqrt(g * (H0 - bmin));
        x = new float[NX];
        z = new float[NZ];
        dx = (xmax - xmin) / NX;
        dz = (zmax - zmin) / NZ;
        for (int i = 0; i < NX; i++) { x[i] = xmin + 0.5f * dx + i * dx; }
        for (int j = 0; j < NZ; j++) { z[j] = zmin + 0.5f * dz + j * dz; }
        drmin = Mathf.Min(dx, dz);

        h = new float[lenArr];
        velx = new float[lenArr];
        velz = new float[lenArr];
        b = new float[lenArr];
        eta = new float[lenArr];
        parbx = new float[lenArr];
        parbz = new float[lenArr];
        fHalf = new Vector3[(NX + 1) * NZ];
        gHalf = new Vector3[NX * (NZ + 1)];

        dispatchX = Mathf.CeilToInt(lenx / (float)THREADS_X);
        dispatchZ = Mathf.CeilToInt(lenz / (float)THREADS_Z);
        dispatch_BDX_X = Mathf.CeilToInt(NGhost / (float)THREADS_X);
        dispatch_BDX_Z = Mathf.CeilToInt(lenz / (float)THREADS_Z);
        dispatch_BDZ_X = Mathf.CeilToInt(lenx / (float)THREADS_X);
        dispatch_BDZ_Z = Mathf.CeilToInt(NGhost / (float)THREADS_Z);

        //SetupSeabedBeach();
        SetupSeabedGaussianRock();
    }

    void SetupSeabedGaussianRock()
    {
        float[] positionX = { 2, 4, 6, 8 };
        float[] positionZ = { 0, -3, 2, 1 };
        float[] hillHight = { 0.5f, 1.0f, 0.8f, 0.7f };
        float bSigma = 0.5f; 
        for (int i = 0; i < lenx; i++)
        {
            for (int j = 0; j < lenz; j++)
            {
                int idnow = i * lenz + j;
                eta[idnow] = 1;
                b[idnow] = bmin;
                parbx[idnow] = 0f;
                parbz[idnow] = 0f;
                if (i >= NGhost && i < NX + NGhost && j >= NGhost && j < NZ + NGhost)
                {
                    float xnow = x[i - NGhost]; 
                    float znow = z[j - NGhost];
                    for (int n=0; n<positionX.Length; n++)
                    {
                        float delx = xnow - positionX[n];
                        float delz = znow - positionZ[n];
                        float r2 = delx * delx + delz * delz;
                        float gaussian = hillHight[n] * Mathf.Exp(-r2 / (2 * bSigma * bSigma));
                        b[idnow] += bmax * gaussian;
                        parbx[idnow] += - bmax * gaussian * delx / (bSigma * bSigma);
                        parbz[idnow] += - bmax * gaussian * delz / (bSigma * bSigma);
                    }
                }
                h[idnow] = Mathf.Max(hmin, eta[idnow] - b[idnow]);
                eta[idnow] = h[idnow] + b[idnow];
                velx[idnow] = 0;
                velz[idnow] = 0;
            }
        }
    }

    void SetupSeabedBeach()
    {
        float kz = bZWaveNumber * 2f * Mathf.PI / (zmax - zmin);
        for (int i = 0; i < lenx; i++)
        {
            for (int j = 0; j < lenz; j++)
            {
                int idnow = i * lenz + j;
                eta[idnow] = 1;
                if (i < NGhost)
                {
                    b[idnow] = bmin;
                    parbx[idnow] = 0f;
                    parbz[idnow] = 0f;

                }
                else if (i >= NX + NGhost)
                {
                    b[idnow] = bmax;
                    parbx[idnow] = 0f;
                    parbz[idnow] = 0f;
                }
                else if (j < NGhost || j >= NZ + NGhost)
                {
                    float slopeX = (bmax / (xmax - xmin)) * (x[i - NGhost] - xmin);
                    b[idnow] = bmin + slopeX * (1f + bAmp);
                    parbx[idnow] = 0f;
                    parbz[idnow] = 0f;
                }
                else
                {
                    float slopeX = (bmax / (xmax - xmin)) * (x[i - NGhost] - xmin);
                    b[idnow] = bmin + slopeX * (1f + bAmp * Mathf.Cos(kz * (z[j - NGhost] - zmin)));
                    parbx[idnow] = (bmax / (xmax - xmin)) * (1f + bAmp * Mathf.Cos(kz * (z[j - NGhost] - zmin)));
                    parbz[idnow] = -slopeX * bAmp * kz * Mathf.Sin(kz * (z[j - NGhost] - zmin));
                }
                h[idnow] = Mathf.Max(hmin, eta[idnow] - b[idnow]);
                eta[idnow] = h[idnow] + b[idnow];
                velx[idnow] = 0;
                velz[idnow] = 0;
            }
        }
    }

    void BuiltMesh()
    {
        Mesh waterMesh = CreateMesh(0);
        MeshFilter mfW = GetComponent<MeshFilter>();
        MeshRenderer mrW = GetComponent<MeshRenderer>();
        mrW.material = WaterMaterial;
        mfW.mesh = waterMesh;

        GameObject seabedGO = new GameObject("Seabed");
        seabedGO.transform.SetParent(transform, false);
        MeshFilter mfB = seabedGO.AddComponent<MeshFilter>();
        MeshRenderer mrB = seabedGO.AddComponent<MeshRenderer>();
        mfB.mesh = CreateMesh(1);
        mrB.material = SeaBedMaterial;
    }

    Mesh CreateMesh(int meshType)
    {
        Mesh mesh = new Mesh();
        mesh.indexFormat = UnityEngine.Rendering.IndexFormat.UInt32;
        Vector3[] vertices = new Vector3[NX * NZ];
        Vector2[] uv = new Vector2[NX * NZ];
        int[] triangles = new int[(NX - 1) * (NZ - 1) *6];

        for (int j = 0; j < NZ; j++) {
            for (int i = 0; i < NX; i++){
                int idnow = i * NZ + j;
                int idMesh = i * NZ + j;
                int idData = (i + NGhost) * lenz + (j + NGhost);
                if (meshType == 0)
                {
                    vertices[idMesh] = new Vector3( xmin + i * dx,
                                                    eta[idData],
                                                    zmin + j * dz );
                }
                else
                {
                    vertices[idMesh] = new Vector3( xmin + i * dx,
                                                    b[idData] + 10*hmin,
                                                    zmin + j * dz );
                }
                    uv[idMesh] = new Vector2((float)i / (NX - 1), (float)j / (NZ - 1));
            }
        }
        int tcount = 0;
        for (int j = 0; j < NZ - 1; j++){
            for (int i = 0; i < NX - 1; i++){
                int idMesh = i * NZ + j;
                triangles[tcount++] = idMesh; 
                triangles[tcount++] = idMesh + 1; 
                triangles[tcount++] = idMesh + NZ; 
                triangles[tcount++] = idMesh + 1; 
                triangles[tcount++] = idMesh + NZ + 1; 
                triangles[tcount++] = idMesh + NZ;
            }
        }
        mesh.vertices = vertices;
        mesh.uv = uv;
        mesh.triangles = triangles;
        mesh.RecalculateNormals();
        mesh.bounds = new Bounds(
            new Vector3((xmin + xmax) / 2f, H0, (zmin + zmax) / 2f),
            new Vector3((xmax - xmin), 10f, (zmax - zmin))
        );
        return mesh;
    }

    void FindKernel()
    {
        kernelUpdateBDX = SWSimShader.FindKernel("kernelUpdateBDX");
        kernelUpdateBDZ = SWSimShader.FindKernel("kernelUpdateBDZ");
        kernelUpdateF = SWSimShader.FindKernel("kernelUpdateF");
        kernelUpdateG = SWSimShader.FindKernel("kernelUpdateG");
        kernelUpdateU = SWSimShader.FindKernel("kernelUpdateU");
    }

    void CreateBuffer()
    {
        hBuffer = new ComputeBuffer(lenArr, sizeof(float));
        velxBuffer = new ComputeBuffer(lenArr, sizeof(float)); 
        velzBuffer = new ComputeBuffer(lenArr, sizeof(float)); 
        bBuffer = new ComputeBuffer(lenArr, sizeof(float)); 
        etaBuffer = new ComputeBuffer(lenArr, sizeof(float));
        parbxBuffer = new ComputeBuffer(lenArr, sizeof(float)); 
        parbzBuffer = new ComputeBuffer(lenArr, sizeof(float));
        fHalfBuffer = new ComputeBuffer((NX + 1) * NZ, sizeof(float) * 3);
        gHalfBuffer = new ComputeBuffer(NX * (NZ + 1), sizeof(float) * 3);
    }

    void bindBuffer()
    {
        SWSimShader.SetBuffer(kernelUpdateBDX, "h", hBuffer);
        SWSimShader.SetBuffer(kernelUpdateBDX, "velx", velxBuffer);
        SWSimShader.SetBuffer(kernelUpdateBDX, "velz", velzBuffer);

        SWSimShader.SetBuffer(kernelUpdateBDZ, "h", hBuffer);
        SWSimShader.SetBuffer(kernelUpdateBDZ, "velx", velxBuffer);
        SWSimShader.SetBuffer(kernelUpdateBDZ, "velz", velzBuffer);

        SWSimShader.SetBuffer(kernelUpdateF, "h", hBuffer);
        SWSimShader.SetBuffer(kernelUpdateF, "velx", velxBuffer);
        SWSimShader.SetBuffer(kernelUpdateF, "velz", velzBuffer);
        SWSimShader.SetBuffer(kernelUpdateF, "fHalf", fHalfBuffer);

        SWSimShader.SetBuffer(kernelUpdateG, "h", hBuffer);
        SWSimShader.SetBuffer(kernelUpdateG, "velx", velxBuffer);
        SWSimShader.SetBuffer(kernelUpdateG, "velz", velzBuffer);
        SWSimShader.SetBuffer(kernelUpdateG, "gHalf", gHalfBuffer);

        SWSimShader.SetBuffer(kernelUpdateU, "h", hBuffer);
        SWSimShader.SetBuffer(kernelUpdateU, "velx", velxBuffer);
        SWSimShader.SetBuffer(kernelUpdateU, "velz", velzBuffer);
        SWSimShader.SetBuffer(kernelUpdateU, "b", bBuffer);
        SWSimShader.SetBuffer(kernelUpdateU, "eta", etaBuffer);
        SWSimShader.SetBuffer(kernelUpdateU, "parbx", parbxBuffer);
        SWSimShader.SetBuffer(kernelUpdateU, "parbz", parbzBuffer);
        SWSimShader.SetBuffer(kernelUpdateU, "fHalf", fHalfBuffer);
        SWSimShader.SetBuffer(kernelUpdateU, "gHalf", gHalfBuffer);
    }

    void SetBufferData()
    {
        hBuffer.SetData(h);
        velxBuffer.SetData(velx);
        velzBuffer.SetData(velz);
        bBuffer.SetData(b);
        etaBuffer.SetData(eta);
        parbxBuffer.SetData(parbx);
        parbzBuffer.SetData(parbz);
    }

    void SetCBuffer()
    {
        // compute shader
        SWSimShader.SetInt("NX", NX); 
        SWSimShader.SetInt("NZ", NZ); 
        SWSimShader.SetInt("NGhost", NGhost); 
        SWSimShader.SetInt("lenx", lenx); 
        SWSimShader.SetInt("lenz", lenz);
        SWSimShader.SetInt("lenArr", lenArr);
        SWSimShader.SetFloat("eta0", H0);
        SWSimShader.SetFloat("g", g); 
        SWSimShader.SetFloat("dx", dx); 
        SWSimShader.SetFloat("dz", dz);
        SWSimShader.SetFloat("bmin", bmin);
        SWSimShader.SetFloat("injectWaveFreq", injectWaveFreq);
        SWSimShader.SetFloat("injectWaveAmp", injectWaveAmp);
        
        //SWSimShader.SetFloat("cmax", cmax);
        SWSimShader.SetFloat("hmin", hmin);
        // render material 
        WaterMaterial.SetInt("NX", NX);
        WaterMaterial.SetInt("NZ", NZ);
        WaterMaterial.SetInt("NGhost", NGhost);
        WaterMaterial.SetInt("lenx", lenx);
        WaterMaterial.SetInt("lenz", lenz);
        WaterMaterial.SetInt("lenArr", lenArr);
        WaterMaterial.SetFloat("dx", dx);
        WaterMaterial.SetFloat("dz", dz);
    }

    void SafeRelease(ref ComputeBuffer buffer)
    {
        if (buffer != null)
        {
            buffer.Release();
            buffer = null;
        }
    }

    void OnDestroy()
    {
        SafeRelease(ref etaBuffer);
        SafeRelease(ref hBuffer);
        SafeRelease(ref velxBuffer);
        SafeRelease(ref velzBuffer);
        SafeRelease(ref bBuffer);
        SafeRelease(ref parbxBuffer);
        SafeRelease(ref parbzBuffer);
        SafeRelease(ref fHalfBuffer);
        SafeRelease(ref gHalfBuffer);
    }
}
