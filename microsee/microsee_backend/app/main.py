"""
app/main.py

FastAPI application entry point.
Registers routers, configures CORS, sets up logging.
"""

from __future__ import annotations

import logging
from contextlib import asynccontextmanager

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

from app.routes import parse

# ── Logging ───────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(name)s  %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


# ── Lifespan ──────────────────────────────────────────────────────────────────
@asynccontextmanager
async def lifespan(app: FastAPI):
    logger.info("MicroSee backend starting up")
    yield
    logger.info("MicroSee backend shutting down")


# ── App ───────────────────────────────────────────────────────────────────────
app = FastAPI(
    title="MicroSee API",
    description=(
        "Backend for the MicroSee microbiome visualisation suite. "
        "Provides QIIME2 file parsing, diversity statistics, and "
        "chart-ready data preparation."
    ),
    version="1.0.0",
    docs_url="/docs",      # Swagger UI
    redoc_url="/redoc",    # ReDoc
    lifespan=lifespan,
)

# ── CORS ──────────────────────────────────────────────────────────────────────
# Allow the HTML file served from any origin (file://, localhost, etc.)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],          # tighten this in production
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ── Routers ───────────────────────────────────────────────────────────────────
app.include_router(parse.router)


# ── Health check ─────────────────────────────────────────────────────────────
@app.get("/health", tags=["health"], summary="Backend health check")
async def health():
    return {"status": "ok", "version": app.version}


# ── Global exception handler ──────────────────────────────────────────────────
@app.exception_handler(Exception)
async def unhandled_exception(request, exc):
    logger.exception(f"Unhandled error on {request.url}: {exc}")
    return JSONResponse(
        status_code=500,
        content={"detail": "Internal server error", "type": type(exc).__name__},
    )
