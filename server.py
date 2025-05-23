from mcp.server.fastmcp import FastMCP
from mcp.server.fastmcp.prompts import base

mcp = FastMCP("scanpy_mcp")

if __name__ == "__main__":
    import uvicorn
    uvicorn.run("main:mcp", host="127.0.0.1", port=8000)
