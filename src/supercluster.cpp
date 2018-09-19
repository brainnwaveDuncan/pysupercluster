/*
    Copyright (c) 2016, Wemap SAS

    Permission to use, copy, modify, and/or distribute this software for any
    purpose with or without fee is hereby granted, provided that the above
    copyright notice and this permission notice appear in all copies.

    THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
    WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
    MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
    ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
    WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
    ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
    OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#include <algorithm>

#include "supercluster.hpp"


Cluster::Cluster(const Point &_point, size_t _numPoints, size_t _id, int _expansionZoom, std::vector<size_t> children)
    : point(_point)
    , numPoints(_numPoints)
    , id(_id)
    , zoom(999)
    , expansionZoom(_expansionZoom)
    , children(children)
{
}


ClusterTree::ClusterTree(const std::vector<Cluster*> &_clusters)
    : clusters(_clusters)
{
    std::vector<Point> points;
    points.reserve(clusters.size());
    for (size_t i = 0; i < clusters.size(); ++i)
        points.push_back(clusters[i]->point);
    kdbush = new kdbush::KDBush<Point>(points);
}


ClusterTree::~ClusterTree()
{
    delete kdbush;
}


SuperCluster::SuperCluster(const std::vector<Point> &points, int _minZoom, int _maxZoom, double _radius, double _extent)
    : minZoom(_minZoom)
    , maxZoom(_maxZoom)
    , radius(_radius)
    , extent(_extent)
{
    trees.resize(maxZoom + 2);

    // prepare initial clusters
    std::vector<Cluster*> clusters;
    clusters.reserve(points.size());
    for (size_t i = 0; i < points.size(); ++i)
        clusters.push_back(new Cluster(points[i], 1, i, -1, {std::get<2>(points[i])}));
    all_clusters = clusters;

    for (int z = maxZoom; z >= minZoom; --z) {
        trees[z + 1] = new ClusterTree(clusters);
        clusters = cluster(clusters, z);
    }

    // index top-level clusters
    trees[minZoom] = new ClusterTree(clusters);
}


SuperCluster::~SuperCluster()
{
    for (size_t i = 0; i < trees.size(); ++i)
        delete trees[i];
    for (size_t i = 0; i < all_clusters.size(); ++i)
        delete all_clusters[i];
}


std::vector<Cluster*> SuperCluster::cluster(const std::vector<Cluster*> &points, int zoom)
{
    std::vector<Cluster*> clusters;
    double radius = this->radius / (this->extent * (1 << zoom));
    ClusterTree *tree = trees[zoom + 1];

    for (size_t i = 0; i < points.size(); ++i) {
        Cluster *p = points[i];
        if (p->zoom <= zoom)
            continue;
        p->zoom = zoom;

        bool foundNeighbors = false;
        size_t numPoints = p->numPoints;
        double wx = std::get<0>(p->point) * numPoints;
        double wy = std::get<1>(p->point) * numPoints;
        std::vector<size_t> children(p->children.begin(), p->children.end());

        tree->kdbush->within(std::get<0>(p->point), std::get<1>(p->point), radius, [&children, &foundNeighbors, &numPoints, tree, &wx, &wy, zoom](const std::vector<Cluster*>::size_type id) {
            Cluster *b = tree->clusters[id];
            if (zoom < b->zoom) {
                foundNeighbors = true;
                b->zoom = zoom;
                wx += std::get<0>(b->point) * b->numPoints;
                wy += std::get<1>(b->point) * b->numPoints;
                numPoints += b->numPoints;
                children.insert(children.end(), b->children.begin(), b->children.end());
            }
        });

        if (foundNeighbors) {
            Cluster *cluster = new Cluster(Point(wx / numPoints, wy / numPoints, 0), numPoints, all_clusters.size(), zoom + 1, children);
            clusters.push_back(cluster);
            all_clusters.push_back(cluster);
        } else {
            clusters.push_back(p);
        }
    }

    return clusters;
}


std::vector<Cluster*> SuperCluster::getClusters(const Point &min_p, const Point &max_p, int zoom) const
{
    const int z = std::max(minZoom, std::min(zoom, maxZoom + 1));
    std::vector<Cluster*> clusters;

    ClusterTree *tree = trees[z];
    tree->kdbush->range(std::get<0>(min_p), std::get<1>(min_p), std::get<0>(max_p), std::get<1>(max_p), [&clusters, &tree](const std::vector<Cluster*>::size_type id) {
        clusters.push_back(tree->clusters[id]);
    });

    return clusters;
}
