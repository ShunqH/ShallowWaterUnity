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

            // TEXTURE2D(_BaseMap);
            // SAMPLER(sampler_BaseMap);

            // CBUFFER_START(UnityPerMaterial)
            //     half4 _BaseColor;
            //     float4 _BaseMap_ST;
            // CBUFFER_END

            Interpolators vert(MeshData IN)
            {
                Interpolators OUT;
                OUT.vertex = TransformObjectToHClip(IN.vertex.xyz);
                float normalStrength = 5.0;
                OUT.normal = IN.normal; 
                OUT.normal = normalize(float3(IN.normal.x * normalStrength, IN.normal.y, IN.normal.z * normalStrength)); 
                OUT.uv = IN.uv; 
                return OUT;
            }

            float4 frag(Interpolators IN) : SV_Target
            {
                float3 n = normalize(IN.normal);
                float3 lightDir = normalize(float3(0.3, 0.5, 0.2));
                float ndl = saturate(dot(n, lightDir));

                float3 sandColor = float3(0.75, 0.70, 0.55);

                float3 color = sandColor * (0.2 + 0.8 * ndl);

                return float4(color, 1.0);
            }
            ENDHLSL
        }
    }
}
