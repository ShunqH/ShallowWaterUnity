Shader "Custom/RenderShader"
{
    Properties
    {
        [MainColor] _BaseColor("Base Color", Color) = (1, 1, 1, 1)
        [MainTexture] _BaseMap("Base Map", 2D) = "white"
    }

    SubShader
    {
        Tags { "RenderType" = "Opaque" "RenderPipeline" = "UniversalPipeline" }

        Pass
        {
            HLSLPROGRAM

            #pragma vertex vert
            #pragma fragment frag

            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Core.hlsl"

            StructuredBuffer<float> eta;

            cbuffer SimulationParams : register(b0)
            {
                int NX, NZ, NGhost;
                int lenx, lenz, lenArr; 
                float dx, dz;
            }

            struct MeshData
            {
                float4 vertex : POSITION;
                float3 normal : NORMAL; 
                float2 uv : TEXCOORD0;
                uint vid : SV_VertexID; 
            };

            struct Interpolators
            {
                float4 vertex : SV_POSITION;
                float3 normal : TEXCOORD0;
                float2 uv : TEXCOORD1; 
            };

            float ComputeEta(uint i, uint j)
            {
                return eta[i*lenz+j]; 
            }

            Interpolators vert(MeshData IN)
            {
                Interpolators OUT;
                uint inow = IN.vid / NZ + NGhost; 
                uint jnow = IN.vid % NZ + NGhost; 

                float H = ComputeEta(inow, jnow); 
                float Hx = ( ComputeEta(inow+1, jnow)  - ComputeEta(inow-1, jnow) ) / ( 2 * dx ); 
                float Hz = ( ComputeEta(inow, jnow+1)  - ComputeEta(inow, jnow-1) ) / ( 2 * dz ); 
                float normalStrength = 20.0;
                float3 n = normalize(float3(-Hx * normalStrength, 1.0, Hz * normalStrength)); 

                IN.vertex.y = H;

                OUT.vertex = TransformObjectToHClip(IN.vertex.xyz);
                OUT.normal = TransformObjectToWorldNormal(n); 
                OUT.uv = IN.uv; 
                return OUT;
            }

            float4 frag(Interpolators IN) : SV_Target
            {
                float3 lightDir = normalize(float3(0.3, 0.5, 0.2));
                float ndl = saturate(dot(normalize(IN.normal), lightDir));
                ndl = pow(ndl, 0.4);
                float3 base = float3(0.0, 0.6, 0.8);
                float3 color = base * (0.15 + 0.85 * ndl);
                return float4(color, 0.8);
            }
            ENDHLSL
        }
    }
}
